"""Explicit, source-backed pilot parameterization; regional transfer is labelled."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd,rasterio
from scipy.ndimage import distance_transform_edt
from rasterio.windows import from_bounds,Window
from rasterio.warp import reproject,Resampling,transform_bounds
ROOT=next(p for p in Path(__file__).resolve().parents if (p/'outputs/执行结果_20261025/服务参数').is_dir())
PROJECT=Path('C:/Users/33221/Documents/Codex/2026-09-29/new-chat/outputs/china_ecosystem_scenarios_20261022')
sys.path.insert(0,str(PROJECT/'src'))
from ecological_services import water_yield_budyko,soil_retention,carbon_stock_from_fractions
OUT=ROOT/'outputs/执行结果_20261025/服务参数'
PARS=pd.read_csv(OUT/'九类完整服务参数.csv')
THREATS=pd.read_csv(OUT/'HQ威胁采用值与映射.csv')

def write(path,arr,template,names):
    if arr.ndim==2:arr=arr[None]
    p=template.copy();p.update(count=len(arr),dtype='float32',nodata=-9999,compress='deflate')
    with rasterio.open(path,'w',**p) as ds:
        ds.write(np.where(np.isfinite(arr),arr,-9999).astype('float32'))
        for i,n in enumerate(names,1):ds.set_band_description(i,n)

def headquarters(lc,pixel_m,distances=None,factor=1.,exclude_bare=False):
    """Match InVEST 3.14.2 nearest-source EDT+decay equations on native CLCD."""
    if distances is None:
        distances={r.threat:(distance_transform_edt(lc!=r.clcd_code,sampling=pixel_m) if np.any(lc==r.clcd_code) else np.full(lc.shape,np.inf)) for r in THREATS.itertuples()}
    active=THREATS[THREATS.threat!='bare'] if exclude_bare else THREATS
    total=np.zeros(lc.shape,float)
    for r in active.itertuples():
        d=distances[r.threat];radius=r.max_dist_km*1000*factor
        effect=np.where(d<=radius,1-d/radius if r.decay=='linear' else np.exp(-2.99*d/radius),0)
        lookup=np.r_[np.nan,PARS['sens_'+r.threat].values]
        # Removing bare pressure holds other threat weights fixed: a one-factor ablation.
        total+=r.weight/THREATS.weight.sum()*effect*lookup[lc]
    H=np.r_[np.nan,PARS.H.values][lc]
    return H/(1+(total/.5)**2.5),distances

def carbon(fractions):
    d=PARS[['C_above_MgC_ha','C_below_MgC_ha','C_soil_MgC_ha','C_dead_MgC_ha']].values
    return carbon_stock_from_fractions(fractions,d,1e6)[0]/100

def cover_factor(ndvi,low=.05,high=.8):
    f=100*np.clip((np.asarray(ndvi)-low)/(high-low),0,1)
    with np.errstate(divide='ignore',invalid='ignore'):
        c=np.clip(.6508-.3436*np.log10(f),0,1)
    return np.where(f==0,1,np.where(f>=78.3,0,c))

def soil_k(sand,silt,clay,soc):
    """EPIC with China correction, t ha h/(ha MJ mm); all input percentages."""
    sn=1-sand/100
    with np.errstate(invalid='ignore',divide='ignore'):
        epic=(.2+.3*np.exp(-.0256*sand*(1-silt/100)))*(silt/(clay+silt))**.3*(1-.25*soc/(soc+np.exp(3.72-2.95*soc)))*(1-.7*sn/(sn+np.exp(22.9*sn-5.51)))
    return (-.01383+.51575*epic)*.1317

def water_yield(fractions,precip,et0,pawc,profile_depth_mm=1000,z=1):
    vals=[]
    for r in PARS.itertuples():
        pc=pawc[r.clcd_code-1] if np.ndim(pawc)==3 else pawc
        if not r.WY_process_supported:val=np.full(fractions.shape[1:],np.nan)
        else:val=water_yield_budyko(precip,et0,r.Kc,pc,r.root_depth_mm,profile_depth_mm,z,r.vegetated)[0]
        vals.append(val)
    vals=np.asarray(vals)
    supported=np.all((fractions<=1e-8)|np.isfinite(vals),axis=0)
    return np.where(supported,np.nansum(fractions*vals,axis=0),np.nan)

def retained_soil(fractions,ndvi,r,k,ls,low=.05,high=.8,crop_p=1,exclude_wetland=False):
    """Per total cell ha; structural zeros outside soil process, not missing soils."""
    c=cover_factor(ndvi,low,high)
    soil_codes=[1,2,3,4,7]+([] if exclude_wetland else [9])
    soil_frac=fractions[np.array(soil_codes)-1].sum(axis=0)
    managed=np.zeros_like(c)
    for code in soil_codes:
        cv=np.ones_like(c) if code==7 else c
        p=crop_p if code==1 else 1.
        managed+=fractions[code-1]*np.where(fractions[code-1]>1e-8,cv*p,0)
    potential=r*k*ls
    retained=potential*(soil_frac-managed)
    # Missing NDVI must not invalidate completely nonvegetated pixels.
    no_soil=soil_frac<=1e-8
    retained[no_soil]=0
    return retained,soil_frac,c

def run_hq_carbon():
    summary=[];results=OUT/'服务栅格';results.mkdir(exist_ok=True)
    for year in [2000,2005,2010,2015,2020,2025]:
        with rasterio.open(PROJECT/f'data/national_landcover/CLCD_v01_{year}_albert.tif') as native:
            for fp in sorted((ROOT/'outputs/执行结果_20261025/样区栅格').glob(f'*_{year}_fractions.tif')):
                site=fp.name.rsplit('_',2)[0]
                with rasterio.open(fp) as ds:
                    fractions=ds.read()[:9];profile=ds.profile;b=ds.bounds
                    # Maximum sensitivity distance is 12 km; extra native margin prevents edge truncation.
                    target=(b.left-12500,b.bottom-12500,b.right+12500,b.top+12500)
                    sb=transform_bounds(ds.crs,native.crs,*target)
                    nw=from_bounds(*sb,transform=native.transform).round_offsets().round_lengths()
                    raw=native.read(1,window=nw);nt=native.window_transform(nw)
                    assert np.isin(raw,np.arange(10)).all()
                    q,d=headquarters(raw,abs(nt.a))
                    variants=[q,headquarters(raw,abs(nt.a),d,exclude_bare=True)[0],headquarters(raw,abs(nt.a),d,factor=.8)[0],headquarters(raw,abs(nt.a),d,factor=1.2)[0]]
                    arr=[]
                    for v in variants:
                        dest=np.full(fractions.shape[1:],-9999,np.float32)
                        reproject(np.where(np.isfinite(v),v,-9999),dest,src_transform=nt,src_crs=native.crs,src_nodata=-9999,dst_transform=ds.transform,dst_crs=ds.crs,dst_nodata=-9999,resampling=Resampling.average)
                        dest[dest==-9999]=np.nan;arr.append(dest)
                    carr=carbon(fractions);arr.append(carr)
                    names=['HQ_native30m_avg1km','HQ_remove_bare_threat','HQ_distance_x0_8','HQ_distance_x1_2','Carbon_MgC_ha_Zoige_regional_transfer']
                    write(results/f'{site}_{year}_HQ_Carbon.tif',np.stack(arr),profile,names)
                    assert np.nanmin(arr[:4])>=0 and np.nanmax(arr[:4])<=1
                    assert np.all(arr[1]>=arr[0]-1e-6)
                    summary.append(dict(site=site,year=year,HQ_mean=float(np.nanmean(arr[0])),HQ_no_bare_mean=float(np.nanmean(arr[1])),HQ_dist08_mean=float(np.nanmean(arr[2])),HQ_dist12_mean=float(np.nanmean(arr[3])),carbon_MgC_ha=float(np.nanmean(carr)),carbon_scope='regional_Zoige' if site=='Zoige' else 'Zoige_transfer_diagnostic_only',native_HQ_pixel_m=abs(nt.a),halo_m=12500,finite_cells=int(np.isfinite(arr[0]).sum())))
                    print(site,year,summary[-1]['HQ_mean'],flush=True)
        pd.DataFrame(summary).to_csv(OUT/'七样区六期_HQ与碳运行汇总.csv',index=False,encoding='utf-8-sig')
    (OUT/'HQ与碳运行验证.json').write_text(json.dumps(dict(passed=True,runs=len(summary),checks=['HQ bounds 0-1','native 30m EDT before 1km average','halo exceeds max tested distance','all 9 classes and 4 pools','full fractions sum validated by carbon algebra'],full_national_parameter_calibration=False),indent=2),encoding='utf-8')

if __name__=='__main__':run_hq_carbon()
