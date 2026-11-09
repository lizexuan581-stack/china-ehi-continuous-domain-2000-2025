"""Explicit current-source four-block accounting; no legacy global imports."""
import ast,hashlib,math,csv
from pathlib import Path
import numpy as np
KEYS=('EV','EO','HQ','Carbon','SR','WY','FP')
GROUPS=('land_direct','meteorology_direct','observed_vegetation','food_map_composite')
class EndpointMismatch(AssertionError):
 def __init__(self,label,rebuilt,expected):
  bad=np.flatnonzero(rebuilt!=expected);self.details=dict(label=label,mismatch_count=len(bad),max_abs_difference=float(np.nanmax(abs(rebuilt.astype(float)-expected.astype(float)))),samples=[dict(local_index=int(i),expected=float(expected[i]),rebuilt=float(rebuilt[i])) for i in bad[:20]]);super().__init__(str(self.details))
def source_helpers(spec):
 ns={'np':np}
 for item in spec['pure_helper_sources']:
  path=Path(item['path']);assert hashlib.sha256(path.read_bytes()).hexdigest()==item['sha256'];tree=ast.parse(path.read_text(encoding='utf-8-sig'));selected=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in item['functions']];assert len(selected)==len(item['functions']);exec(compile(ast.Module(body=selected,type_ignores=[]),str(path),'exec'),ns)
 return ns['water_yield_budyko'],ns['cover_factor']
def parameters(path):
 with Path(path).open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
 assert [int(x['clcd_code']) for x in rows]==list(range(1,10));return rows
def normalize(k,a,anchors):
 a=np.asarray(a,'float64')
 if k in anchors:a=(a-anchors[k]['low'])/(anchors[k]['high']-anchors[k]['low'])
 return np.clip(a,0,1)
def score(ev,eo,hq,carbon,sr,wy,food):return np.cbrt(ev*eo*hq)*(carbon+sr+wy+food)/4
def water(lc,p,et,static,pars,helper):
 # Annual production reads all physical arrays as float32; helper casts to f64.
 p=np.asarray(p,'float32');et=np.asarray(et,'float32');out=np.full(lc.shape,np.nan,'float64')
 for row in pars:
  code=int(row['clcd_code']);sel=lc==code
  if not sel.any() or int(row['WY_process_supported'])==0:continue
  pa=static['pawet'] if code==9 else static['pa'];v=helper(p[sel],et[sel],float(row['Kc']),pa[sel],float(row['root_depth_mm']),static['bed'][sel],1,int(row['vegetated']))[0];assert np.all(~np.isfinite(v)|((v>=-1e-6)&(v<=p[sel]+1e-5)));out[sel]=v
 return out.astype('float32')
def soil(lc,ndvi,static,cover):
 # Preserve original float32 operation order, including empirical cover dtype.
 ndvi=np.asarray(ndvi,'float32');r,k,ls=[np.asarray(static[x],'float32') for x in ('r','k','ls')];cv=cover(ndvi);raw=r*k*ls*(1-cv);out=np.zeros(lc.shape,'float64');veg=np.isin(lc,[1,2,3,4,9]);out[veg]=raw[veg];return out.astype('float32')
def subset_shapley(values):
 assert values.shape[0]==16 and np.isfinite(values).all();out=np.zeros((4,values.shape[1]),'float64')
 for j in range(4):
  for state in range(16):
   if not state&(1<<j):
    n=state.bit_count();out[j]+=math.factorial(n)*math.factorial(3-n)/24*(values[state|(1<<j)]-values[state])
 return out
def evaluate(endpoints,static,anchors,pars,helpers,check_physical_endpoints=True):
 water_helper,cover=helpers;n=len(endpoints[0]['EV']);wy={};sr={};errors={}
 for d in endpoints:
  assert all(np.isfinite(d[k]).all() for k in KEYS) and np.isfinite(d['P']).all() and np.isfinite(d['ET0']).all()
 for l in (0,1):
  for pi in (0,1):
   for ei in (0,1):
    if l==pi==ei:
     expected=np.asarray(endpoints[l]['WY'],'float32')
     if check_physical_endpoints:
      rebuilt=water(endpoints[l]['land'],endpoints[pi]['P'],endpoints[ei]['ET0'],static,pars,water_helper);errors['WY_'+str(l)]=float(np.max(abs(rebuilt.astype(float)-expected.astype(float))))
      if not np.array_equal(rebuilt,expected):raise EndpointMismatch('WY_'+str(l),rebuilt,expected)
     wy[l,pi,ei]=expected
    else:wy[l,pi,ei]=water(endpoints[l]['land'],endpoints[pi]['P'],endpoints[ei]['ET0'],static,pars,water_helper)
  for v in (0,1):
   if l==v:
    expected=np.asarray(endpoints[l]['SR'],'float32')
    if check_physical_endpoints:
     rebuilt=soil(endpoints[l]['land'],endpoints[v]['EV'],static,cover);errors['SR_'+str(l)]=float(np.max(abs(rebuilt.astype(float)-expected.astype(float))))
     if not np.array_equal(rebuilt,expected):raise EndpointMismatch('SR_'+str(l),rebuilt,expected)
    sr[l,v]=expected
   else:sr[l,v]=soil(endpoints[l]['land'],endpoints[v]['EV'],static,cover)
 wy_norm={k:normalize('WY',v,anchors) for k,v in wy.items()};sr_norm={k:normalize('SR',v,anchors) for k,v in sr.items()};norm=[{k:normalize(k,d[k],anchors) for k in KEYS} for d in endpoints]
 values=np.empty((16,n),'float64');hybrid={}
 for state in range(16):
  l=state&1;c=(state>>1)&1;v=(state>>2)&1;f=(state>>3)&1
  def evaluate_weather(pi,ei):return score(norm[v]['EV'],norm[l]['EO'],norm[l]['HQ'],norm[l]['Carbon'],sr_norm[l,v],wy_norm[l,pi,ei],norm[f]['FP'])
  values[state]=evaluate_weather(c,c)
  if not c:hybrid[state,1,0]=evaluate_weather(1,0);hybrid[state,0,1]=evaluate_weather(0,1)
 valid=np.isfinite(values).all(axis=0)
 for h in hybrid.values():valid&=np.isfinite(h)
 assert valid.all(),'Mixed-state invalidity on fixed domain; denominator must not shrink'
 phi=subset_shapley(values);precip=np.zeros(n);demand=np.zeros(n)
 for state in range(16):
  if state&2:continue
  sz=state.bit_count();weight=math.factorial(sz)*math.factorial(3-sz)/24
  precip+=weight*.5*((hybrid[state,1,0]-values[state])+(values[state|2]-hybrid[state,0,1]))
  demand+=weight*.5*((hybrid[state,0,1]-values[state])+(values[state|2]-hybrid[state,1,0]))
 assert np.max(abs(phi.sum(0)-(values[15]-values[0])))<1e-12 and np.max(abs(precip+demand-phi[1]))<1e-12
 return dict(values=values,phi=phi,precipitation_path=precip,ET0_path=demand,wy=wy,sr=sr,physical_endpoint_errors=errors,hybrid=hybrid)
