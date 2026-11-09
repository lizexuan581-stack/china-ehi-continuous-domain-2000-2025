"""Complete-component national EHI; no service dropped or weight reassigned."""
from national_physical_common import *
from contextlib import ExitStack
from rasterio.windows import Window
import pandas as pd,datetime,sys
YEARS=[2000,2005,2010,2015,2020,2025]
KEYS=['EV','EO','HQ','Carbon','SR','WY','FP'];NORMAL=['EV','Carbon','SR','WY','FP']
DEST=D.parent/'national_ehi';DEST.mkdir(exist_ok=True)
VARIANTS={'main':{},'loose_NDVI':{'EV':3,'SR':2,'FP':4},'WY_Z5':{'WY':2},'WY_Z10':{'WY':3},'HQ_remove_bare':{'HQ':2},'SR_original_slope':{'SR':3},'SR_annual_R':{'SR':4},'FP_area':{'FP':1},'FP_observed_conditional':{'FP':3}}
def paths(y):
    return {'EV':D.parent/f'fp_national_cache/EV_{y}_QA01_DETAILED_1km.tif','FP':D.parent/f'fp_national_cache/FP_{y}_national_mean_weight_model.tif',**{k:COMP/f'{k}_{y}_1km.tif' for k in ['EO','HQ','Carbon','SR','WY']}}
def raw(src,w,inside,variant='main'):
    result={};v=VARIANTS[variant]
    for k,s in src.items():
        a=s.read(v.get(k,2 if k=='FP' else 1),window=w,masked=True).filled(np.nan).astype('float64')
        if k=='FP':a/=100
        a[~inside]=np.nan;result[k]=a
    return result
def integrate(a,anc):
    scored={};flags=np.zeros(a['EV'].shape,'uint16')
    for i,k in enumerate(KEYS):
        s=(a[k]-anc[k]['low'])/(anc[k]['high']-anc[k]['low']) if k in NORMAL else a[k]
        flags[np.isfinite(s)&((s<0)|(s>1))]|=1<<i;scored[k]=np.clip(s,0,1)
    complete=np.logical_and.reduce([np.isfinite(s) for s in scored.values()]);phi=np.cbrt(scored['EV']*scored['EO']*scored['HQ'])
    es=(scored['Carbon']+scored['SR']+scored['WY']+scored['HQ']+scored['FP'])/5
    es4=(scored['Carbon']+scored['SR']+scored['WY']+scored['FP'])/4
    ehi=np.where(complete,np.sqrt(phi*es),np.nan);ehi4=np.where(complete,np.sqrt(phi*es4),np.nan)
    assert np.all(~np.isfinite(ehi)|((ehi>=0)&(ehi<=1)))
    return ehi,phi,es,ehi4,scored,flags,complete
def main():
    weight=json.loads((ROOT/'work/service_parameters/es_weighting.json').read_text());assert weight['weights']==[.2]*5
    worker='--years' in sys.argv
    runyears=[int(y) for y in sys.argv[sys.argv.index('--years')+1].split(',')] if worker else YEARS
    p=template();m=mask();anc={k:dict(low=np.inf,high=-np.inf,reference_valid_cells=0) for k in NORMAL}
    if worker:anc=json.loads((OUT/'全国EHI固定历史归一化锚点.json').read_text(encoding='utf-8'))
    for y in ([] if worker else YEARS[:-1]):
        for k in NORMAL:
            with rasterio.open(paths(y)[k]) as s:
                for offset in range(0,p['height'],512):
                    h=min(512,p['height']-offset);w=Window(0,offset,p['width'],h);a=s.read(2 if k=='FP' else 1,window=w,masked=True).filled(np.nan);good=np.isfinite(a)&m[offset:offset+h]
                    if k=='FP':a=a/100
                    if good.any():anc[k]['low']=min(anc[k]['low'],float(a[good].min()));anc[k]['high']=max(anc[k]['high'],float(a[good].max()));anc[k]['reference_valid_cells']+=int(good.sum())
        print('ANCHORS',y,flush=True)
    assert all(v['high']>v['low'] for v in anc.values())
    if not worker:(OUT/'全国EHI固定历史归一化锚点.json').write_text(json.dumps(anc,ensure_ascii=False,indent=2),encoding='utf-8')
    allstats=[];sources=[]
    with rasterio.open(ROOT/'work/Province_ID_GADM41_1km.tif') as ids:
      for y in runyears:
        with ExitStack() as stack:
            src={k:stack.enter_context(rasterio.open(fn)) for k,fn in paths(y).items()};expected=(p['height'],p['width'],p['transform'],p['crs'])
            assert all((s.height,s.width,s.transform,s.crs)==expected for s in src.values())
            op=p.copy();op.update(dtype='float32',nodata=-9999,compress='deflate',zlevel=1,tiled=True,blockxsize=256,blockysize=256,BIGTIFF='IF_SAFER')
            dst={v:stack.enter_context(rasterio.open(DEST/f'EHI_{y}_{v}.tif','w',count=14 if v=='main' else 1,**{k:x for k,x in op.items() if k!='count'})) for v in VARIANTS}
            stats={v:np.zeros((100,5),float) for v in VARIANTS};missing={k:0 for k in KEYS};valid_any=0
            for offset in range(0,p['height'],256):
                h=min(256,p['height']-offset);w=Window(0,offset,p['width'],h);inside=m[offset:offset+h];pid=ids.read(1,window=w)
                base=raw(src,w,inside,'main')
                for v in VARIANTS:
                    a=base.copy()
                    for k,b in VARIANTS[v].items():
                        value=src[k].read(b,window=w,masked=True).filled(np.nan).astype('float64')
                        if k=='FP':value/=100
                        value[~inside]=np.nan;a[k]=value
                    e,phi,es,e4,s,flags,complete=integrate(a,anc)
                    code=pid[complete];data=[np.ones(len(code)),e[complete],e[complete]**2,(flags[complete]>0).astype(float),(e[complete]==0).astype(float)]
                    for i,vals in enumerate(data):stats[v][:,i]+=np.bincount(code,weights=vals,minlength=100)
                    arrays=[e]
                    if v=='main':
                        fm=src['FP'].read(6,window=w,masked=True).filled(np.nan);arrays=[e,phi,es,e4]+[s[k] for k in KEYS]+[flags.astype(float),fm,complete.astype(float)]
                        for k in KEYS:missing[k]+=int((inside&~np.isfinite(a[k])).sum())
                    for b,arr in enumerate(arrays,1):dst[v].write(np.where(inside&np.isfinite(arr),arr,-9999).astype('float32'),b,window=w)
                if offset%1280==0:print('EHI ROWS',y,offset,p['height'],flush=True)
            for v,s in stats.items():
                for code in range(100):
                    if s[code,0]>0:allstats.append(dict(year=y,variant=v,province_code=code,valid_cells=int(s[code,0]),EHI_sum=s[code,1],EHI_mean=s[code,1]/s[code,0],clipped_cells=int(s[code,3]),zero_EHI_cells=int(s[code,4])))
            names=['EHI_ES5','PHI_VOR','ES5_equal_weights','EHI_ES4_no_HQ_in_services']+[k+'_score' for k in KEYS]+['out_of_reference_bitmask','FP_weight_modeled_mask','complete_component_mask']
            for b,n in enumerate(names,1):dst['main'].set_band_description(b,n)
            for v in VARIANTS:
                if v!='main':dst[v].set_band_description(1,'EHI_ES5_'+v)
                dst[v].update_tags(scope='31 mainland provincial analysis domain; finite complete-component support only',weights='0.2 Carbon SR WY HQ FP',reference='all finite primary component cells 2000-2020; 2025 held out',variant=v)
            sources.append(dict(year=y,paths={k:str(f) for k,f in paths(y).items()},missing_component_cells=missing))
        print('EHI YEAR DONE',y,flush=True)
    resultout=OUT/'分批运行记录' if worker else OUT;resultout.mkdir(exist_ok=True);suffix='_'+str(runyears[0]) if worker else ''
    pd.DataFrame(allstats).to_csv(resultout/('全国EHI各省逐期可用像元汇总'+suffix+'.csv'),index=False,encoding='utf-8-sig')
    config=dict(complete=True,years=YEARS,primary_files=6,variant_count=len(VARIANTS),formula='PHI=cuberoot(EV*EO*HQ); ES=(C+SR+WY+HQ+FP)/5; EHI=sqrt(PHI*ES)',service_weighting=weight,reference_years=YEARS[:-1],heldout=2025,normalization='national fixed historical minmax for EV,C,SR,WY,FP; EO theoretical score; HQ native 0-1; scores clipped with flags',analysis_mask_cells=int(m.sum()),missing_rule='all seven components required; no partial ES mean or missing-as-zero',implementation='1km national categorical land grid; 5km EO; HydroSHEDS terrain; provincial-conserving FP model with flags; this is an auditable national implementation, not exact original Zhu data replication',sources=sources,variants=VARIANTS,timestamp=datetime.datetime.now().astimezone().isoformat())
    config['years']=runyears;config['primary_files']=len(runyears);config['scope']='worker partition' if worker else 'all six dates'
    (resultout/('全国完整EHI运行记录'+suffix+'.json')).write_text(json.dumps(config,ensure_ascii=False,indent=2),encoding='utf-8')
    print('NATIONAL EHI COMPLETE',flush=True)
if __name__=='__main__':main()

