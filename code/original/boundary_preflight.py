from pathlib import Path
import json,csv,hashlib
import numpy as np,rasterio
from rasterio.features import rasterize
from rasterio.warp import transform_geom
from rasterio.enums import MergeAlg
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'outputs/执行结果_20261025'
PROJECT=Path('C:/Users/33221/Documents/Codex/2026-09-29/new-chat/outputs/china_ecosystem_scenarios_20261022')
NAMES={'Beijing':('北京',11),'Tianjin':('天津',12),'Hebei':('河北',13),'Shanxi':('山西',14),'NeiMongol':('内蒙古',15),'Liaoning':('辽宁',21),'Jilin':('吉林',22),'Heilongjiang':('黑龙江',23),'Shanghai':('上海',31),'Jiangsu':('江苏',32),'Zhejiang':('浙江',33),'Anhui':('安徽',34),'Fujian':('福建',35),'Jiangxi':('江西',36),'Shandong':('山东',37),'Henan':('河南',41),'Hubei':('湖北',42),'Hunan':('湖南',43),'Guangdong':('广东',44),'Guangxi':('广西',45),'Hainan':('海南',46),'Chongqing':('重庆',50),'Sichuan':('四川',51),'Guizhou':('贵州',52),'Yunnan':('云南',53),'Xizang':('西藏',54),'Shaanxi':('陕西',61),'Gansu':('甘肃',62),'Qinghai':('青海',63),'NingxiaHui':('宁夏',64),'XinjiangUygur':('新疆',65)}
def main():
    source=ROOT/'sources/gadm41_CHN_1.json';j=json.loads(source.read_text(encoding='utf-8'))
    template=PROJECT/'data/national_landcover/CLCD_2025_ScenarioGrid_1km.tif'
    with rasterio.open(template) as ds:
        shapes=[(transform_geom('EPSG:4326',ds.crs,f['geometry']),NAMES[f['properties']['NAME_1']][1]) for f in j['features'] if f['properties']['NAME_1'] in NAMES]
        ids=rasterize(shapes,out_shape=ds.shape,transform=ds.transform,fill=0,dtype='uint8',all_touched=False)
        overlaps=rasterize([(g,1) for g,c in shapes],out_shape=ds.shape,transform=ds.transform,fill=0,dtype='uint8',merge_alg=MergeAlg.add)
        reverse=rasterize(list(reversed(shapes)),out_shape=ds.shape,transform=ds.transform,fill=0,dtype='uint8')
        conflict=ids!=reverse
        ids[conflict]=0
        profile=ds.profile.copy();profile.update(count=1,dtype='uint8',nodata=0,compress='deflate')
        # Derived administrative raster stays in work because GADM restricts redistribution.
        with rasterio.open(ROOT/'Province_ID_GADM41_1km.tif','w',**profile) as dst:dst.write(ids,1)
        grid=dict(width=ds.width,height=ds.height,transform=list(ds.transform)[:6],crs_wkt=ds.crs.to_wkt(),pixel_area_m2=abs(ds.transform.a*ds.transform.e))
    rows=[];yr=[]
    for year in [2000,2005,2010,2015,2020,2025]:
        with rasterio.open(PROJECT/f'data/national_landcover/CLCD_{year}_ScenarioGrid_1km.tif') as ds:lc=ds.read(1)
        good=(lc>0)&(ids>0)
        yr.append(dict(year=year,clcd_valid_cells=int((lc>0).sum()),boundary_intersection_cells=int(good.sum()),clcd_cells_outside_boundary=int(((lc>0)&(ids==0)).sum()),boundary_cells_without_clcd=int(((ids>0)&(lc==0)).sum())))
        for en,(zh,code) in NAMES.items():
            mask=ids==code
            rows.append(dict(year=year,province=zh,province_code=code,gadm_name=en,boundary_cells=int(mask.sum()),valid_lc_cells=int((mask&(lc>0)).sum()),mode_cropland_cells=int((mask&(lc==1)).sum())))
    with (OUT/'省域边界与土地覆盖核对.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    report=dict(source='GADM4.1 CHN ADM1',source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),source_features=len(j['features']),
        selected_geometry_count=len(shapes),selected_province_names=len(NAMES),repeated_names='Multiple features assigned same province ID; no geometry discarded within selected names',
        exclusion=['HongKong','Macau'],grid=grid,center_rule='pixel center; all_touched=False',overlapping_center_cells=int((overlaps>1).sum()),
        conflicting_province_center_cells=int(conflict.sum()),conflict_rule='exclude conflicting province-center assignments; same-province multipart overlap retained',
        by_year=yr,source_file_retained_locally=True,boundary_status='adopted analysis boundary v1; CLCD intersection explicit',
        map_boundary_note='Analysis boundary and statistical mask; publication cartography prepared separately',raster_path=str(ROOT/'Province_ID_GADM41_1km.tif'))
    (OUT/'边界与网格预检.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['source_features','selected_geometry_count','selected_province_names','overlapping_center_cells','by_year']},ensure_ascii=False))
if __name__=='__main__':main()
