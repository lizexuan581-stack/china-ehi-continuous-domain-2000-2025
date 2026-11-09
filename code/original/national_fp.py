"""National, blockwise crop-specific FP inputs and explicit allocation comparators.

Observed-support NDVI allocation is conditional on observed cropland. Missing
cropland is accounted for, never silently presented as a full-support model.
"""
from pathlib import Path
from contextlib import ExitStack
import sys,json,argparse,time
import numpy as np,pandas as pd,rasterio
from rasterio.windows import Window,from_bounds
from rasterio.warp import transform_bounds,reproject,Resampling
from compose_pilots import manifest,masks,restore_four,PROJECT,OriginalTileRepair,DOYS,BASEOUT,equivalent_grid

ROOT=Path(__file__).resolve().parents[2]
CACHE=Path('D:/Research_Data/China_EHI_20261022/fp_national_cache');CACHE.mkdir(exist_ok=True)
OUT=BASEOUT/'历史植被';IDS=ROOT/'work/Province_ID_GADM41_1km.tif'

def record(path,obj):path.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
def blocks(ds,n=512):
    for r in range(0,ds.height,n):
        for c in range(0,ds.width,n):yield Window(c,r,min(n,ds.width-c),min(n,ds.height-r))
def profile(ds,count):
    p=ds.profile.copy();p.update(driver='GTiff',count=count,dtype='float32',nodata=-9999,compress='deflate',tiled=True,blockxsize=256,blockysize=256,BIGTIFF='IF_SAFER');return p

def native_composite(year):
    dest=CACHE/f'NDVI_{year}_QA01_DETAILED_native.tif';done=dest.with_suffix('.json')
    if dest.exists() and done.exists() and json.loads(done.read_text())['complete']:return dest
    resolved=manifest(year)
    if resolved is None:raise RuntimeError('Require verified complete source year')
    meta,idx=resolved
    with rasterio.Env(GDAL_CACHEMAX=256*1024*1024),ExitStack() as stack:
        src={k:stack.enter_context(rasterio.open(v)) for k,v in idx.items()};base=src['NDVI',49]
        assert all(equivalent_grid(x,base) for x in src.values())
        repair=OriginalTileRepair(PROJECT,base,stack) if year==2025 else None
        ds=stack.enter_context(rasterio.open(dest,'w',**profile(base,5)))
        for j,w in enumerate(blocks(base)):
            sh=(int(w.height),int(w.width));s=np.zeros(sh,np.float64);n=np.zeros(sh,np.uint8);g=np.zeros((4,*sh),bool)
            for k,d in enumerate(DOYS):
                nd=src['NDVI',d].read(1,window=w);rel=src['pixel_reliability',d].read(1,window=w)
                day=src['composite_day_of_the_year',d].read(1,window=w);qa=src['VI_Quality',d].read(1,window=w)
                if repair:nd,rel,day,qa=restore_four(repair,nd,rel,day,qa,d,w)
                good=masks(nd,rel,day,qa,d,year)[0]['QA01_DETAILED']
                s+=np.where(good,nd*.0001,0);n+=good;g[k//5]|=good
            gc=g.sum(axis=0);mean=np.divide(s,n,out=np.full(sh,-9999.),where=n>0)
            ds.write(np.stack([mean,n,gc,(n>=10)&(gc>=4),(n>=5)&(gc>=3)]).astype('float32'),window=w)
            if j%20==0:print(year,'native block',j,flush=True)
        for i,name in enumerate(['NDVI_observed_mean','valid_composite_count','observed_time_groups','qualified_n10g4','qualified_n5g3'],1):ds.set_band_description(i,name)
        rr=repair.report() if repair else {'used':False}
    record(done,{'complete':True,'year':year,'rule':'QA01_DETAILED','source_restoration':rr,'source_complete':True,'temporal_interpolation':False})
    return dest

def crop_components(year,native_path):
    from crop_parallel import run
    return run(year,native_path,CACHE,IDS,PROJECT)

def crop_components_serial_reference(year,native_path):
    dest=CACHE/f'Crop_{year}_30m_NDVI_components_1km.tif';done=dest.with_suffix('.json')
    if dest.exists() and done.exists() and json.loads(done.read_text())['complete']:return dest
    tile_cache=CACHE/f'crop_blocks_{year}';tile_cache.mkdir(exist_ok=True)
    with rasterio.Env(GDAL_CACHEMAX=256*1024*1024),ExitStack() as stack:
        grid=stack.enter_context(rasterio.open(IDS));lc=stack.enter_context(rasterio.open(PROJECT/f'data/national_landcover/CLCD_v01_{year}_albert.tif'))
        nv=stack.enter_context(rasterio.open(native_path));ds=stack.enter_context(rasterio.open(dest,'w',**profile(grid,8)))
        assert abs(abs(grid.transform.a*grid.transform.e)-1e6)<1e-6
        processed=0
        for j,w in enumerate(blocks(grid,128)):
            sh=(int(w.height),int(w.width));a=np.zeros((8,*sh),np.float32);pid=grid.read(1,window=w)
            if not (pid>0).any():ds.write(a,window=w);continue
            checkpoint=tile_cache/f'{j:05d}.npz'
            if checkpoint.exists():
                with np.load(checkpoint) as saved:a=saved['components']
                assert a.shape==(8,*sh)
                ds.write(a,window=w);processed+=1;continue
            tf=grid.window_transform(w);bounds=grid.window_bounds(w)
            cb=transform_bounds(grid.crs,lc.crs,*bounds);cw=from_bounds(*cb,transform=lc.transform).round_offsets().round_lengths();cw=Window(cw.col_off-3,cw.row_off-3,cw.width+6,cw.height+6)
            raw=lc.read(1,window=cw,boundless=True,fill_value=0);crop=raw==1;ct=lc.window_transform(cw)
            finekw=dict(src_transform=ct,src_crs=lc.crs,dst_transform=tf,dst_crs=grid.crs,resampling=Resampling.average,num_threads=4)
            reproject(crop.astype('float32'),a[0],**finekw);reproject((raw>0).astype('float32'),a[7],**finekw)
            if crop.any():
                nb=transform_bounds(grid.crs,nv.crs,*bounds);nw=from_bounds(*nb,transform=nv.transform).round_offsets().round_lengths();nw=Window(nw.col_off-4,nw.row_off-4,nw.width+8,nw.height+8)
                values=nv.read([1,4,5],window=nw,boundless=True,fill_value=-9999);nt=nv.window_transform(nw)
                for k in range(2):
                    valid=(values[k+1]==1)&(values[0]>=-.2)&(values[0]<=1)
                    fine=np.full(crop.shape,-9999,np.float32)
                    reproject(np.where(valid,values[0],-9999).astype('float32'),fine,src_transform=nt,src_crs=nv.crs,src_nodata=-9999,dst_transform=ct,dst_crs=lc.crs,dst_nodata=-9999,resampling=Resampling.nearest,num_threads=4)
                    known=crop&(fine>=0)&(fine<=1);neg=crop&(fine<0)&(fine>=-.2)
                    reproject(known.astype('float32'),a[1+k*3],**finekw)
                    reproject(np.where(known,fine,0).astype('float32'),a[2+k*3],**finekw)
                    reproject(neg.astype('float32'),a[3+k*3],**finekw)
            assert np.all(a[1]<=a[0]+1e-6) and np.all(a[4]<=a[0]+1e-6)
            temp=checkpoint.with_suffix('.tmp.npz');np.savez_compressed(temp,components=a);temp.replace(checkpoint)
            ds.write(a,window=w);processed+=1
            if processed%10==0:print(year,'crop block',processed,'grid block',j,flush=True)
        names=['cropland_area_km2','observed_cropland_area_n10g4_km2','observed_crop_NDVI_weight_n10g4_km2','negative_NDVI_crop_area_n10g4_km2','observed_cropland_area_n5g3_km2','observed_crop_NDVI_weight_n5g3_km2','negative_NDVI_crop_area_n5g3_km2','native_landcover_coverage_fraction']
        for i,name in enumerate(names,1):ds.set_band_description(i,name)
    record(done,{'complete':True,'year':year,'native_NDVI_cache':str(native_path),'crop_class':1,'cropland_resolution_m':30,'allocation_grid_m':1000,'negative_NDVI_policy':'exclude from qualified nonnegative weighting support; area recorded separately','bands':names})
    return dest

def allocate(year,component_path):
    grain=pd.read_csv(BASEOUT/'六期省级粮食_统一字段.csv');grain=grain[grain.year==year]
    with rasterio.open(component_path) as ds,rasterio.open(IDS) as ids:
        sums=np.zeros((100,8),np.float64)
        for _,w in ds.block_windows(1):
            pid=ids.read(1,window=w);v=ds.read(window=w)
            for i in range(8):sums[:,i]+=np.bincount(pid.ravel(),weights=v[i].ravel(),minlength=100)[:100]
        rows=[];factors=np.zeros((100,3),np.float64)
        for r in grain.itertuples():
            code=r.province_code;area=sums[code,0];assert area>0
            factors[code,0]=r.production_tonnes/area
            for k,(oi,wi,ni) in enumerate([(1,2,3),(4,5,6)]):
                denom=sums[code,wi];assert denom>0
                factors[code,k+1]=r.production_tonnes/denom
                rows.append(dict(year=year,province=r.province,province_code=code,temporal_support=['n10g4','n5g3'][k],source_production_tonnes=r.production_tonnes,cropland_area_km2=area,qualified_nonnegative_crop_area_km2=sums[code,oi],qualified_crop_area_pct=100*sums[code,oi]/area,negative_NDVI_crop_area_km2=sums[code,ni],observed_NDVI_area_weight_km2=denom,allocation_scope='conditional_on_observed_crop_support',full_crop_support=bool(sums[code,oi]/area>=.999999)))
        dest=CACHE/f'FP_{year}_national_allocation_comparators.tif';realized=np.zeros((100,3),np.float64)
        with rasterio.open(dest,'w',**profile(ds,3)) as out:
            for _,w in ds.block_windows(1):
                pid=ids.read(1,window=w);v=ds.read(window=w);a=np.stack([v[0]*factors[pid,0],v[2]*factors[pid,1],v[5]*factors[pid,2]])
                for k in range(3):realized[:,k]+=np.bincount(pid.ravel(),weights=a[k].ravel(),minlength=100)[:100]
                # Partial-support weight is an observed contribution, not estimated missing-crop production.
                a[:,pid==0]=-9999;out.write(a.astype('float32'),window=w)
            for i,name in enumerate(['FP_true_30m_crop_area_only_tonnes','FP_NDVI_n10g4_observed_support_conditional_tonnes','FP_NDVI_n5g3_observed_support_conditional_tonnes'],1):out.set_band_description(i,name)
        checks=[]
        # Verify persisted float32 output, including quantization, by province.
        persisted=np.zeros((100,3),np.float64)
        with rasterio.open(dest) as final:
            for _,w in final.block_windows(1):
                pid=ids.read(1,window=w);a=final.read(window=w);a[:,pid==0]=0
                for k in range(3):persisted[:,k]+=np.bincount(pid.ravel(),weights=a[k].ravel(),minlength=100)[:100]
            for r in grain.itertuples():
                for k,method in enumerate(['true_crop_area_only','NDVI_n10g4_observed_support','NDVI_n5g3_observed_support']):
                    total=persisted[r.province_code,k];err=(total-r.production_tonnes)/r.production_tonnes;assert abs(err)<1e-6
                    checks.append(dict(year=year,province=r.province,province_code=r.province_code,method=method,source_tonnes=r.production_tonnes,persisted_allocated_tonnes=total,relative_error=err))
            for p in sorted((BASEOUT/'样区栅格').glob(f'*_{year}_fractions.tif')):
                site=p.name.split(f'_{year}_')[0]
                with rasterio.open(p) as target:
                    w=from_bounds(*target.bounds,transform=final.transform).round_offsets().round_lengths();a=final.read(window=w)
                    # Include crop coverage so conditional allocation is not confused with complete FP.
                    c=ds.read([1,2,5],window=w);a=np.concatenate([a,c])
                    with rasterio.open(OUT/f'{site}_{year}_FP_national_denominator_comparators.tif','w',**profile(target,6)) as trial:
                        trial.write(a);names=list(final.descriptions)+['cropland_area_km2','observed_crop_area_n10g4_km2','observed_crop_area_n5g3_km2']
                        for i,name in enumerate(names,1):trial.set_band_description(i,name)
        pd.DataFrame(rows).to_csv(OUT/f'FP_{year}_省级NDVI分母与覆盖.csv',index=False,encoding='utf-8-sig')
        pd.DataFrame(checks).to_csv(OUT/f'FP_{year}_三方案省级守恒.csv',index=False,encoding='utf-8-sig')
        record(OUT/f'FP_{year}_全国分配运行记录.json',{'complete':True,'year':year,'national_component_cache':str(component_path),'allocation_cache':str(dest),'provinces':len(grain),'pilot_slices':7,'max_relative_conservation_error':max(abs(r['relative_error']) for r in checks),'native_30m_crop_area_used':True,'NDVI_model_status':'conditional observed-support comparator; inspect provincial missing crop share before primary use','national_source_revision_harmonized':False})
    print(year,'FP national allocation finished',flush=True)

def allocate_mean_model(year,component_path):
    """Model only missing allocation weights; do not alter any observed NDVI."""
    grain=pd.read_csv(BASEOUT/'六期省级粮食_统一字段.csv');grain=grain[grain.year==year]
    with rasterio.open(component_path) as ds,rasterio.open(IDS) as ids:
        sums=np.zeros((100,8),np.float64)
        for _,w in ds.block_windows(1):
            pid=ids.read(1,window=w);v=ds.read(window=w)
            for i in range(8):sums[:,i]+=np.bincount(pid.ravel(),weights=v[i].ravel(),minlength=100)[:100]
        rows=[];factors=np.zeros((100,5),np.float64);means=np.zeros((100,2),np.float64);fallback=np.zeros((100,2),bool)
        for r in grain.itertuples():
            code=r.province_code;area=sums[code,0];assert area>0
            factors[code,0]=r.production_tonnes/area
            for k,(oi,wi,ni) in enumerate([(1,2,3),(4,5,6)]):
                obs=sums[code,oi]+sums[code,ni];missing=max(area-obs,0);sobs=sums[code,wi]
                mean=sobs/obs if obs>0 else 0.;means[code,k]=mean
                modeled=sobs+missing*mean
                fallback[code,k]=(modeled<=0)
                factors[code,1+2*k]=r.production_tonnes/(modeled if modeled>0 else area)
                factors[code,2+2*k]=r.production_tonnes/(sobs if sobs>0 else area)
                rows.append(dict(year=year,province=r.province,province_code=code,temporal_support=['n10g4','n5g3'][k],source_production_tonnes=r.production_tonnes,cropland_area_km2=area,observed_crop_area_km2=obs,missing_crop_area_km2=missing,missing_crop_pct=100*missing/area,negative_NDVI_crop_area_km2=sums[code,ni],observed_nonnegative_NDVI_area_weight_km2=sobs,province_observed_area_weighted_mean_NDVI=mean,modeled_missing_area_weight_km2=missing*mean,modeled_total_weight_km2=modeled,theoretical_weight_lower=sobs,theoretical_weight_upper=sobs+missing,area_fallback=bool(fallback[code,k]),fallback_reason='no_positive_observed_weight' if fallback[code,k] else '',negative_NDVI_weight_rule='zero weight; observed area retained',missing_weight_rule='same-province observed-crop area-weighted mean; not an observed NDVI value'))
        dest=CACHE/f'FP_{year}_national_mean_weight_model.tif'
        names=['FP_true_30m_crop_area_only_tonnes','FP_mean_weight_model_n10g4_tonnes','FP_observed_support_conditional_n10g4_tonnes','FP_mean_weight_model_n5g3_tonnes','FP_observed_support_conditional_n5g3_tonnes','missing_weight_model_mask_n10g4','missing_weight_model_mask_n5g3','missing_crop_fraction_n10g4','missing_crop_fraction_n5g3','negative_NDVI_crop_area_n10g4_km2','negative_NDVI_crop_area_n5g3_km2']
        diff=np.zeros((100,4),np.float64)
        with rasterio.open(dest,'w',**profile(ds,11)) as out:
            for _,w in ds.block_windows(1):
                pid=ids.read(1,window=w);v=ds.read(window=w);a=np.zeros((11,*pid.shape),np.float64);a[0]=v[0]*factors[pid,0]
                for k,(oi,wi,ni) in enumerate([(1,2,3),(4,5,6)]):
                    missing=np.maximum(v[0].astype('float64')-v[oi]-v[ni],0)
                    modelw=v[wi]+missing*means[pid,k]
                    modelw=np.where(fallback[pid,k],v[0],modelw)
                    conditionw=np.where(fallback[pid,k],v[0],v[wi])
                    a[1+2*k]=modelw*factors[pid,1+2*k];a[2+2*k]=conditionw*factors[pid,2+2*k]
                    a[5+k]=(missing>1e-6)&(pid>0);a[7+k]=np.divide(missing,v[0],out=np.zeros_like(missing),where=v[0]>0);a[9+k]=v[ni]
                    for j,reference in enumerate([a[0],a[2+2*k]]):diff[:,2*k+j]+=np.bincount(pid.ravel(),weights=np.abs(a[1+2*k]-reference).ravel(),minlength=100)[:100]
                a[:,pid==0]=-9999;out.write(a.astype('float32'),window=w)
            for i,name in enumerate(names,1):out.set_band_description(i,name)
            out.update_tags(quantity='estimated grain supply in tonnes per 1km cell, not measured pixel yield',missing_weight_model='same-province observed crop area-weighted mean; raw NDVI untouched',negative_weight_rule='max(observed NDVI,0); negative crop area retained and reported')
        persisted=np.zeros((100,5),np.float64);checks=[]
        with rasterio.open(dest) as final:
            for _,w in final.block_windows(1):
                pid=ids.read(1,window=w);a=final.read([1,2,3,4,5],window=w);a[:,pid==0]=0
                assert np.all(a>=0)
                for k in range(5):persisted[:,k]+=np.bincount(pid.ravel(),weights=a[k].ravel(),minlength=100)[:100]
            for r in grain.itertuples():
                for k in range(5):
                    total=persisted[r.province_code,k];err=(total-r.production_tonnes)/r.production_tonnes;assert abs(err)<1e-6
                    checks.append(dict(year=year,province=r.province,province_code=r.province_code,method=names[k],source_tonnes=r.production_tonnes,persisted_allocated_tonnes=total,relative_error=err))
            for p in sorted((BASEOUT/'样区栅格').glob(f'*_{year}_fractions.tif')):
                site=p.name.split(f'_{year}_')[0]
                with rasterio.open(p) as target:
                    w=from_bounds(*target.bounds,transform=final.transform).round_offsets().round_lengths();a=final.read(window=w)
                    with rasterio.open(OUT/f'{site}_{year}_FP_national_mean_weight_model.tif','w',**profile(target,11)) as trial:
                        trial.write(a)
                        for i,name in enumerate(names,1):trial.set_band_description(i,name)
        for row in rows:
            code=row['province_code'];k=0 if row['temporal_support']=='n10g4' else 1
            row['mean_model_vs_area_L1_tonnes']=diff[code,2*k]
            row['mean_model_vs_conditional_L1_tonnes']=diff[code,2*k+1]
            row['mean_model_vs_area_redistributed_share']=.5*diff[code,2*k]/row['source_production_tonnes']
            row['mean_model_vs_conditional_redistributed_share']=.5*diff[code,2*k+1]/row['source_production_tonnes']
        pd.DataFrame(rows).to_csv(OUT/f'FP_{year}_省级均值模型与敏感性.csv',index=False,encoding='utf-8-sig')
        pd.DataFrame(checks).to_csv(OUT/f'FP_{year}_均值模型三类方案守恒.csv',index=False,encoding='utf-8-sig')
        record(OUT/f'FP_{year}_均值权重建模运行记录.json',{'complete':True,'year':year,'component_cache':str(component_path),'allocation_cache':str(dest),'bands':names,'provinces':len(grain),'pilot_slices':7,'negative_NDVI_rule':'nonnegative allocation weight max(NDVI,0); raw NDVI unchanged','missing_weight_model':'same-province observed cropland area-weighted mean; area fallback if no positive observed weight','model_value_interpretation':'estimated food provision; not remotely sensed or measured pixel production','max_relative_conservation_error':max(abs(r['relative_error']) for r in checks),'grain_revision_harmonized':False})
    print(year,'FP mean-weight model finished',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--years',type=int,nargs='+',default=[2000,2025]);p.add_argument('--allocation-only',action='store_true');args=p.parse_args()
    for year in args.years:
        # A separately scheduled native-only job may own this year's TIFF.
        # Wait for its completion receipt before opening or replacing the file.
        claim=CACHE/f'native_precompute_{year}.lock'
        if claim.exists():
            deadline=time.monotonic()+3600
            receipt=CACHE/f'NDVI_{year}_QA01_DETAILED_native.json'
            while not (receipt.exists() and json.loads(receipt.read_text()).get('complete')):
                if time.monotonic()>deadline:raise RuntimeError(f'Native precompute {year} did not finish; inspect its log before resuming')
                time.sleep(10)
            claim.unlink(missing_ok=True)
        if args.allocation_only:c=CACHE/f'Crop_{year}_30m_NDVI_components_1km.tif'
        else:n=native_composite(year);c=crop_components(year,n)
        allocate_mean_model(year,c)
