from pathlib import Path
import sys,json,time,csv,hashlib,os
import numpy as np
import rasterio
ROOT=Path(__file__).resolve().parents[2];WORK=Path(__file__).resolve().parent
OUT=ROOT/'results_run';RES=OUT/'results';CACHE=ROOT/'work/cache'
RES.mkdir(parents=True,exist_ok=True);CACHE.mkdir(parents=True,exist_ok=True)
load=lambda p:json.loads(Path(p).read_text(encoding='utf-8-sig'))
m=load(ROOT/'metadata/manifest.json')
source_map=load(os.environ['EHI_SOURCE_MAP']) if os.environ.get('EHI_SOURCE_MAP') else {}
def absolute(obj):
 if isinstance(obj,dict):return {k:absolute(v) for k,v in obj.items()}
 if isinstance(obj,list):return [absolute(v) for v in obj]
 if isinstance(obj,str) and obj.startswith(('data/','metadata/','code/')):return source_map.get(obj,str(ROOT/obj))
 return obj
m=absolute(m);src=m['source_binding'];DEST=Path(src['common26']).parent
cfg={'years':{},'FP_registry':{}}
for y,v in m['years'].items():
 cfg['years'][y]={k:{'path':p} for k,p in v['components'].items()}
 cfg['FP_registry'][y]={'output':{'path':v['components']['FP']}}
anc=load(src['anchors'])['anchors'];YEARS=np.arange(2000,2026);N=4072765
sys.path.insert(0,str(ROOT/'code/original'));import component_core as core;import process_core as process
water_helper,_=process.source_helpers(src);pars=process.parameters(src['parameter_csv'])
def read(path,band=1):
 with rasterio.open(path) as d:return d.read(band,masked=True).astype('float32').filled(np.nan)
def rawpath(y,k):return cfg['FP_registry'][str(y)]['output']['path'] if k=='FP' else cfg['years'][str(y)][k]['path']
def landpath(y):return m['years'][str(y)]['land']
def climatepath(y):return DEST.parent/f'annual_climate_1km/POWER_{y}_strict_1km.tif'
def raw(y,k,mask):return read(rawpath(y,k),2 if k=='FP' else 1)[mask].astype('float64')/(100 if k=='FP' else 1)
def norm(k,x,bounds=anc):return np.clip((x-bounds[k]['low'])/(bounds[k]['high']-bounds[k]['low']),0,1) if k in bounds else np.clip(x,0,1)
def static(mask):return {k:read(v['path'],v['band'])[mask] for k,v in src['static'].items() if k in ['bed','pa','pawet']}
def water(lc,p,e,st,z):
 out=np.full(lc.shape,np.nan,'float32')
 for row in pars:
  code=int(row['clcd_code']);sel=lc==code
  if not sel.any() or not int(row['WY_process_supported']):continue
  pa=st['pawet'] if code==9 else st['pa']
  out[sel]=water_helper(np.asarray(p[sel],'float32'),np.asarray(e[sel],'float32'),float(row['Kc']),pa[sel],float(row['root_depth_mm']),st['bed'][sel],z,int(row['vegetated']))[0].astype('float32')
 return out
def csvout(name,rows):
 fields=list(dict.fromkeys(k for row in rows for k in row))
 with (RES/(name+'.csv')).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
  for row in rows:w.writerow({k:format(v,'.4g') if isinstance(v,(float,np.floating)) else v for k,v in row.items()})
 (CACHE/(name+'_full_precision.json')).write_text(json.dumps(rows,ensure_ascii=False,indent=2,default=lambda x:x.item()),encoding='utf8')
def record(task,t0,**extra):
 value=dict(task=task,seconds=time.monotonic()-t0,**extra)
 (OUT/(task+'_audit.json')).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=lambda x:x.item()),encoding='utf8');print(json.dumps(value,ensure_ascii=False),flush=True)
def sen(y):return float(np.median([(y[j]-y[i])/(j-i) for i in range(len(y)) for j in range(i+1,len(y))]))
def mismatch(a,b,ndvia,ndvib):
 return {'raw':int((core.classify(ndvib-ndvia,b-a,1e-6,1e-6)==2).sum()),'screened':int((core.classify(ndvib-ndvia,b-a,.005,.001)==2).sum())}
def plot_setup():
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 plt.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
 return plt
def savefig(fig,name):
 fig.savefig(RES/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(CACHE/(name+'.pdf'),bbox_inches='tight')
