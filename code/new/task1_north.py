"""Overlay published map line; explicitly retain its cartographic provenance.
The source includes a western extension. Closing segments beyond endpoints are
reported separately; this is an operational geographic partition, not a climate law.
"""
from common import *
import pandas as pd
from pyproj import Transformer
from shapely.geometry import Polygon
from shapely import contains_xy
t0=time.monotonic();mask=np.load(CACHE/'main_mask.npy');ref=np.load(CACHE/'reference_mask.npy')
d=pd.read_stata(WORK/'sources/chinaprov40_line_coord.dta');v=d[d.group==4].dropna();assert set(v._ID)=={38}
crs='+proj=lcc +lat_1=30 +lat_2=62 +lat_0=0 +lon_0=105 +x_0=0 +y_0=0 +ellps=krass +units=m +no_defs'
t=Transformer.from_crs(crs,'EPSG:4326',always_xy=True);lon,lat=t.transform(v._X.values,v._Y.values)
# Source travels east to west. Close via northern exterior of the raster mask.
poly=Polygon(list(zip(lon,lat))+[(60,float(lat[-1])),(60,60),(150,60),(150,float(lat[0]))]);assert poly.is_valid
with rasterio.open(src['common26']) as ds:
 trans=ds.transform;xy=Transformer.from_crs(ds.crs,'EPSG:4326',always_xy=True);lineproj=Transformer.from_crs(crs,ds.crs,always_xy=True)
 north=np.zeros(mask.shape,bool);outside_endpoint=np.zeros(mask.shape,bool)
 for r in range(0,mask.shape[0],256):
  rr,cc=np.where(ref[r:r+256]);xx=trans.c+(cc+.5)*trans.a;yy=trans.f+(rr+r+.5)*trans.e;x,y=xy.transform(xx,yy)
  north[rr+r,cc]=contains_xy(poly,x,y);outside_endpoint[rr+r,cc]=(x<min(lon))|(x>max(lon))
 xline,yline=lineproj.transform(v._X.values,v._Y.values);col=(np.array(xline)-trans.c)/trans.a;row=(np.array(yline)-trans.f)/trans.e
rows=load(CACHE/'task1_domain_composition_full_precision.json')
for label,selected in [('north',north),('south',~north)]:
 r=ref&selected;c=mask&selected;rows.append(dict(group='Qinling_Huaihe_operational_partition',region=label,reference_cells=int(r.sum()),main_cells=int(c.sum()),coverage_pct=100*c.sum()/r.sum(),share_of_main_pct=100*c.sum()/N,status='complete_cartographic_boundary',endpoint_extension_main_cells=int((c&outside_endpoint).sum())))
csvout('task1_domain_composition',rows)
np.save(CACHE/'north_mask.npy',north)
geo={'type':'Feature','properties':{'source':'https://github.com/Wsf921109/stata-prov-map','source_layer':'group=4, ID=38','original_crs':crs,'interpretation':'Cartographic Qinling-Huaihe line including published western extension; full-domain partition closes horizontally beyond endpoints.'},'geometry':{'type':'LineString','coordinates':list(map(list,zip(lon,lat)))}}
(OUT/'Qinling_Huaihe_line_source.geojson').write_text(json.dumps(geo,indent=2),encoding='utf8')
plt=plot_setup();from matplotlib.colors import ListedColormap;from matplotlib.patches import Patch;from matplotlib.lines import Line2D
zone=np.load(CACHE/'aridity_zones.npy');fig,ax=plt.subplots(figsize=(7.2,6.7));show=np.where(mask,zone,np.where(ref,5,0));ax.imshow(show,cmap=ListedColormap(['white','#d8b365','#f6e8c3','#80cdc1','#35978f','#dedede']),vmin=0,vmax=5,interpolation='nearest');ax.plot(col,row,color='#363636',lw=.9);ax.set_axis_off()
handles=[Patch(color=c,label=l) for c,l in zip(['#d8b365','#f6e8c3','#80cdc1','#35978f','#dedede'],['Arid (AI < 0.2)','Semi-arid (0.2–0.5)','Dry subhumid (0.5–0.65)','Humid (AI ≥ 0.65)','Outside continuous observation domain'])]+[Line2D([],[],color='#363636',label='Published Qinling–Huaihe map line')]
ax.legend(handles=handles,loc='lower left',fontsize=8);ax.set_title('Observation support and geographic partitions',loc='left');savefig(fig,'task1_domain_composition');plt.close(fig)
audit=load(OUT/'task1_audit.json');audit.update(north_south_status='complete with explicit cartographic boundary',north_main_cells=int((north&mask).sum()),north_main_pct=100*(north&mask).sum()/N,endpoint_extension_main_cells=int((mask&outside_endpoint).sum()),line_source_sha256=hashlib.sha256((WORK/'sources/chinaprov40_line_coord.dta').read_bytes()).hexdigest(),north_overlay_seconds=time.monotonic()-t0)
(OUT/'task1_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf8');print(audit,flush=True)
