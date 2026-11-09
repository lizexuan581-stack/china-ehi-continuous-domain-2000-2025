"""Regenerate the restricted GADM administrative raster from its official source.
GADM permits noncommercial academic use but restricts redistribution. The source
vector and administrative ID raster are therefore not shipped in this release.
"""
from common import *
import argparse,urllib.request,zipfile
from rasterio.features import rasterize
from rasterio.warp import transform_geom
parser=argparse.ArgumentParser();parser.add_argument('--gadm-zip',type=Path);args=parser.parse_args()
meta=load(ROOT/'metadata/province_regeneration.json');mapping=load(ROOT/'metadata/province_mapping.json')
p=args.gadm_zip or CACHE/'gadm41_CHN_1.json.zip'
if not p.exists():urllib.request.urlretrieve('https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_CHN_1.json.zip',p)
assert hashlib.sha256(p.read_bytes()).hexdigest()==meta['zip_sha256'],'GADM source version differs'
with zipfile.ZipFile(p) as z:
 name=next(n for n in z.namelist() if n.endswith('gadm41_CHN_1.json'));geo=json.loads(z.read(name))
with rasterio.open(absolute(meta['grid_template'])) as ds:
 shapes=[(transform_geom('EPSG:4326',ds.crs,f['geometry']),mapping[f['properties']['NAME_1']][1]) for f in geo['features'] if f['properties']['NAME_1'] in mapping]
 ids=rasterize(shapes,out_shape=ds.shape,transform=ds.transform,fill=0,dtype='uint8',all_touched=False)
 reverse=rasterize(list(reversed(shapes)),out_shape=ds.shape,transform=ds.transform,fill=0,dtype='uint8');ids[ids!=reverse]=0
 assert hashlib.sha256(ids.tobytes()).hexdigest()==meta['array_sha256'],'Rebuilt province pixels differ from original'
 profile=ds.profile.copy();profile.update(count=1,dtype='uint8',nodata=0,compress='deflate')
 target=Path(src['province']);target.parent.mkdir(parents=True,exist_ok=True)
 with rasterio.open(target,'w',**profile) as out:out.write(ids,1)
print('Province raster regenerated; exact original pixel-array checksum passed.')
