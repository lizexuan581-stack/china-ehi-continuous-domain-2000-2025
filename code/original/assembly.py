"""Isolated annual FP and ES4 product assembly. Nothing runs on import.

Reuses the algebra of historical_composition/national_fp.py, national_ehi.py
and es4_revision/historical_es4.py without importing their six-date side effects.
Prepare/preflight are metadata-only. Production needs separate hash-bound review.
"""
from pathlib import Path
from contextlib import ExitStack
import argparse,csv,hashlib,json,os,shutil,time
import numpy as np
import rasterio
from rasterio.windows import Window

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
DATA=Path('D:/Research_Data/China_EHI_20261022/major_revision_20261028')
ASSEMBLY=HERE.parent/'annual_assembly_readiness_20261030'
PHYSICAL=HERE.parent/'annual_physical_20261030'
YEARS=list(range(2000,2026));REFERENCE=[2000,2005,2010,2015,2020]
CODES=[11,12,13,14,15,21,22,23,31,32,33,34,35,36,37,41,42,43,44,45,46,50,51,52,53,54,61,62,63,64,65]
KEYS=['EV','EO','HQ','Carbon','SR','WY','FP'];NORMAL=['EV','Carbon','SR','WY','FP']
SUPPORTS=['n10g4','n5g3']
EV_NAMES=['NDVI_n10g4','qualified_native_area_fraction_n10g4','NDVI_n5g3','qualified_native_area_fraction_n5g3']
CROP_NAMES=['cropland_area_km2','observed_cropland_area_n10g4_km2','observed_crop_NDVI_weight_n10g4_km2','negative_NDVI_crop_area_n10g4_km2','observed_cropland_area_n5g3_km2','observed_crop_NDVI_weight_n5g3_km2','negative_NDVI_crop_area_n5g3_km2','native_landcover_coverage_fraction']
FP_NAMES=['FP_area_tonnes','FP_mean_n10g4_tonnes','FP_observed_n10g4_tonnes','FP_mean_n5g3_tonnes','FP_observed_n5g3_tonnes','missing_model_n10g4','missing_model_n5g3','missing_fraction_n10g4','missing_fraction_n5g3','negative_area_n10g4_km2','negative_area_n5g3_km2']
EHI_NAMES=['EHI_ES4_product','PHI','ES4','EHI_ES4_sqrt','EHI_ES5_sqrt']+[k+'_score' for k in KEYS]+['out_of_reference_bitmask','complete_component_mask']
PHYSICAL_NAMES={'EO':['EO_eight_metrics_5km_8neighbor'],'Carbon':['Carbon_MgC_ha','Biomass_AG_BG_combined','DOM','SOC_fixed_2000_2005'],'HQ':['HQ_fixed_historical_max_degradation_K','HQ_remove_bare'],'WY':['WY_mm_yr_Z1','WY_mm_yr_Z5','WY_mm_yr_Z10','WY_mm_yr_Z20'],'SR':['SR_t_ha_yr_n10g4_Rref','SR_t_ha_yr_n5g3_Rref','SR_original_slope_n10g4_Rref','SR_annual_R_n10g4','soil_process_mask']}
FORMULA='PHI=cuberoot(EV*EO*HQ); ES4=(Carbon+SR+WY+FP)/4; EHI=PHI*ES4'

class GateError(RuntimeError):pass
def require(ok,message):
 if not ok:raise GateError(message)
def load(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def dump(p,obj):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(f'.{os.getpid()}.{time.time_ns()}.tmp.json');tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8');atomic_replace(tmp,p)
def atomic_replace(source,target,attempts=20):
 for attempt in range(attempts):
  try:os.replace(source,target);return
  except PermissionError:
   if attempt==attempts-1:raise
   time.sleep(.1)
def stage_event(dest,status,**details):
 dest=Path(dest);p=dest/'status.json';state=load(p) if p.exists() else {}
 state.update(pid=os.getpid(),status=status,time=time.strftime('%Y-%m-%dT%H:%M:%S'),**details);dump(p,state)
 line=json.dumps(dict(time=state['time'],status=status,**details),ensure_ascii=False)+'\n'
 for attempt in range(20):
  try:
   with (dest/'stages.jsonl').open('a',encoding='utf-8') as f:f.write(line)
   break
  except PermissionError:
   if attempt==19:raise
   time.sleep(.1)
def fp(p):
 p=Path(p);s=p.stat();return dict(path=str(p.resolve()),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def check_fp(v):
 a=fp(v['path']);require(all(a[k]==v[k] for k in a),'Fingerprint changed: '+v['path'])
def walk_bindings(obj):
 if isinstance(obj,dict):
  if {'path','bytes','mtime_ns'}<=set(obj):check_fp(obj)
  for v in obj.values():walk_bindings(v)
 elif isinstance(obj,list):
  for v in obj:walk_bindings(v)
def windows(ds):
 for r in range(0,ds.height,256):
  for c in range(0,ds.width,256):yield Window(c,r,min(256,ds.width-c),min(256,ds.height-r))
def read(ds,band,w):return ds.read(band,window=w,masked=True).astype('float64').filled(np.nan)
def profile(ds,n):
 p=ds.profile.copy();p.update(driver='GTiff',count=n,dtype='float32',nodata=-9999,compress='deflate',predictor=3,tiled=True,blockxsize=256,blockysize=256,BIGTIFF='IF_SAFER');return p

def validate_crop(v):
 require(v.shape[0]==8 and np.isfinite(v).all() and (v>=0).all() and (v<=1.00001).all(),'Invalid crop schema/values')
 for oi,wi,ni in ((1,2,3),(4,5,6)):
  require((v[wi]<=v[oi]+1e-6).all() and (v[oi]+v[ni]<=v[0]+1e-6).all(),'Crop area/weight inequality')

def accumulate_crop(sums,missing_sums,pid,v):
 """Exact per-cell missing area, avoiding global clipping vs local clipping drift."""
 validate_crop(v);require(np.isin(pid,[0]+CODES).all(),'Unexpected province code')
 for i in range(8):sums[:,i]+=np.bincount(pid.ravel(),weights=v[i].ravel(),minlength=100)
 for j,(oi,ni) in enumerate(((1,3),(4,6))):missing_sums[:,j]+=np.bincount(pid.ravel(),weights=np.maximum(v[0]-v[oi]-v[ni],0).ravel(),minlength=100)

def allocation_parameters(sums,missing_sums,budgets):
 factors=np.zeros((100,5));means=np.zeros((100,2));fallback=np.zeros((100,2),bool);rows=[]
 for code,tonnes in budgets.items():
  area=sums[code,0];require(np.isfinite(tonnes) and tonnes>=0 and area>0,'Invalid budget/crop area '+str(code))
  factors[code,0]=tonnes/area
  for j,(oi,wi,ni) in enumerate(((1,2,3),(4,5,6))):
   observed=sums[code,oi]+sums[code,ni];weight=sums[code,wi];missing=missing_sums[code,j]
   mean=weight/observed if observed>0 else 0.;means[code,j]=mean;modeled=weight+missing*mean
   fallback[code,j]=modeled<=0
   factors[code,1+2*j]=tonnes/(modeled if modeled>0 else area)
   factors[code,2+2*j]=tonnes/(weight if weight>0 else area)
   rows.append(dict(province_code=code,support=SUPPORTS[j],production_tonnes=tonnes,crop_area_km2=area,observed_area_km2=observed,missing_area_km2=missing,negative_area_km2=sums[code,ni],mean_weight=mean,modeled_missing_weight=missing*mean,area_fallback=bool(fallback[code,j])))
 return factors,means,fallback,rows

def allocate_block(pid,v,parameters):
 factors,means,fallback,_=parameters;a=np.zeros((11,*pid.shape),float);a[0]=v[0]*factors[pid,0]
 for j,(oi,wi,ni) in enumerate(((1,2,3),(4,5,6))):
  missing=np.maximum(v[0]-v[oi]-v[ni],0);modeled=np.where(fallback[pid,j],v[0],v[wi]+missing*means[pid,j]);conditional=np.where(fallback[pid,j],v[0],v[wi])
  a[1+2*j]=modeled*factors[pid,1+2*j];a[2+2*j]=conditional*factors[pid,2+2*j]
  a[5+j]=(missing>1e-6)&(pid>0);a[7+j]=np.divide(missing,v[0],out=np.zeros_like(missing),where=v[0]>0);a[9+j]=v[ni]
 a[:,pid==0]=np.nan;return a

def empty_anchors():return {k:dict(low=float('inf'),high=-float('inf'),reference_valid_cells=0) for k in NORMAL}
def update_anchors(anc,raw,inside,year):
 require(year in REFERENCE,'Held-out/nonreference year cannot fit anchors')
 for k in NORMAL:
  good=inside&np.isfinite(raw[k]);v=raw[k][good]
  if v.size:anc[k]['low']=min(anc[k]['low'],float(v.min()));anc[k]['high']=max(anc[k]['high'],float(v.max()));anc[k]['reference_valid_cells']+=int(v.size)
def score(raw,anc):
 scored={};flags=np.zeros(raw['EV'].shape,'uint16')
 for i,k in enumerate(KEYS):
  if k in NORMAL:
   require(np.isfinite(anc[k]['low']) and anc[k]['high']>anc[k]['low'],'Degenerate anchors '+k)
   x=(raw[k]-anc[k]['low'])/(anc[k]['high']-anc[k]['low'])
  else:x=raw[k]
  flags[np.isfinite(x)&((x<0)|(x>1))]|=1<<i;scored[k]=np.clip(x,0,1)
 complete=np.logical_and.reduce([np.isfinite(x) for x in scored.values()]);phi=np.cbrt(scored['EV']*scored['EO']*scored['HQ']);es4=sum(scored[k] for k in ('Carbon','SR','WY','FP'))/4;es5=(4*es4+scored['HQ'])/5
 science=[phi*es4,phi,es4,np.sqrt(phi*es4),np.sqrt(phi*es5)]+[scored[k] for k in KEYS]
 return np.stack([np.where(complete,v,np.nan) for v in science]+[flags.astype(float),complete.astype(float)])

def proof(path,checks):
 p=Path(path);return dict(path=str(p.resolve()),sha256=sha(p) if p.exists() else None,checks=checks)
def spec(path,proofs):
 p=Path(path);d=dict(path=str(p.resolve()),fingerprint=fp(p) if p.exists() else None,proofs=proofs)
 if p.exists():
  with rasterio.open(p) as s:d['schema']=list(s.descriptions)
 return d

def prepare(path):
 """Build a reviewable candidate manifest; no rule or budget is promoted."""
 physical=load(PHYSICAL/'annual_physical_queue_status.json');years={}
 for y in YEARS:
  ent={}
  for stage in ('land','WY','SR'):
   q=physical['completed'][f'{y}_{stage}'];pp=[proof(q['support_and_value_audit'],['passed','exact_expected_mask_match'])] if q.get('reused_legacy') else [proof(q['receipt'],['complete','full_read_passed']),proof(q['independent_audit'],['passed','full_read_passed','source_bindings_unchanged']),proof(q['support_audit'],['passed','exact_expected_mask_match']),proof(Path(q['receipt']).with_name('source_binding.json'),[])]
   for k,v in q['outputs'].items():
    if k!='Degradation':ent[k]=spec(v['path'],pp)
  ev=DATA/f'annual_1km_candidates/{y}_ev/ev_{y}_QA01_DETAILED_1km.tif'
  evproof=ev.parent/'complete.json'
  if y==2000:ev=DATA/'annual_1km_candidate_audit/EV_2000_QA01_DETAILED_candidate.tif';evproof=ASSEMBLY/'ev_2000_full_audit.json'
  ent['EV']=spec(ev,[proof(evproof,['passed','complete'] if y==2000 else ['complete','full_read_equality_passed'])])
  crop=DATA/f'annual_crop_candidate_parallel_v1/{y}/Crop_{y}_QA01_DETAILED_candidates_1km.tif'
  if y in REFERENCE+[2025]:
   crop=DATA.parent/f'fp_national_cache/Crop_{y}_30m_NDVI_components_1km.tif'
   pp=[proof(ASSEMBLY/f'legacy_crop_{y}_readonly_audit.json',['complete_readonly_checks','native_full_grid_equivalence','full_crop_checkpoint_equality','values_area_support_constraints_passed']),proof(ASSEMBLY/f'legacy_crop_{y}_kernel_equivalence.json',['complete'])]
  else:pp=[proof(crop.parent/'complete.json',['complete','full_read_equality_passed','area_weight_inequalities_passed']),proof(crop.parent/'independent_audit.json',['complete','full_output_checkpoint_equality','checkpoint_checksum_passed','source_and_code_binding_passed','support_subset_passed','area_weight_constraints_passed']),proof(crop.parent/'source_binding.json',[])]
  ent['crop']=spec(crop,pp)
  # Legacy existing files alone are insufficient: both audits must be complete.
  years[str(y)]=ent
 ids=ROOT/'work/Province_ID_GADM41_1km.tif';land=load(PHYSICAL/'annual_physical_queue_status.json')['completed']['2020_SR']['receipt'];land2020=load(land)['source_binding']['land2020']['path']
 cfg=dict(version='annual_assembly_v1',years=years,grid=spec(ids,[]),land2020=spec(land2020,[]),reference_years=REFERENCE,normalization='fixed_historical_componentwise_minmax',formula=FORMULA,fp_method='same_province_observed_crop_mean_negative_zero_area_fallback',qa_lock=str(ASSEMBLY/'ANNUAL_QA_SUPPORT_POLICY.json'),grain_lock=str(HERE/'GRAIN_FREEZE.json'),output_root=str(DATA/'annual_ehi_assembly_v1'),collection_proof=proof(PHYSICAL/'collection_independent_receipts_audit_20261030.json',['passed']),production_authorization=str(HERE/'PRODUCTION_AUTHORIZATION.json'),authoring_status='candidate_only_not_authorized')
 require(not Path(path).exists(),'Refusing overwrite candidate manifest');dump(path,cfg)

def grain_gate(path,required_years=None):
 required_years=YEARS if required_years is None else list(required_years)
 lock=load(path);require(lock.get('frozen') is True,'Grain version freeze absent')
 require(lock.get('version_id') and lock.get('version_definition') and isinstance(lock.get('unresolved_limitations'),list),'Specific version definition and limitations must be recorded')
 require(lock.get('known_material_production_conflicts_resolved') is True,'Known material production conflicts unresolved')
 require(not lock.get('unresolved_material_version_leads',[]),'Unresolved material production version leads; restrict a pilot lock to unaffected years')
 require(set(required_years)<=set(lock.get('years',[]))<=set(YEARS) and lock.get('province_codes')==CODES,'Grain lock must cover requested years x 31 provinces')
 require(lock.get('units')=='tonnes' and lock.get('allow_interpolation') is False and lock.get('allow_national_rescaling') is False,'Grain units/prohibited substitution')
 if lock.get('scientific_lock'):require(sha(lock['scientific_lock'])==lock.get('scientific_lock_sha256'),'Grain scientific version record changed')
 p=Path(lock['csv']);require(sha(p)==lock['csv_sha256'],'Grain CSV version changed')
 rows=list(csv.DictReader(p.open(encoding='utf-8-sig',newline='')))
 budgets={};seen=set();source_hashes={}
 allowed_status={'CNKI_user_confirmed','official_historical_table','official_report_published_precision','official_contemporaneous_publication'}
 for row in rows:
  y=int(row['year']);code=int(row['province_code']);key=(y,code);require(key not in seen and y in YEARS and code in CODES,'Duplicate/unexpected grain key');seen.add(key)
  if y not in required_years:continue
  t=float(row['production_tonnes']);require(np.isfinite(t) and t>=0,'Invalid grain tonnes')
  require(row.get('source_status') in allowed_status and (row.get('source_edition') or (row.get('version') and row.get('source_identity'))) and len(row['source_sha256'])==64,'Unqualified source status or edition '+str(key))
  require(row.get('known_material_production_conflict','false').lower() not in ('true','1','yes'),'Known material conflict in row '+str(key))
  require(row.get('source_url') or row['source_status']=='CNKI_user_confirmed','Official source URL missing '+str(key))
  evidence=Path(row['source_file']);require(evidence.is_file(),'Missing grain evidence '+str(key))
  if str(evidence) not in source_hashes:source_hashes[str(evidence)]=sha(evidence)
  require(source_hashes[str(evidence)]==row['source_sha256'],'Grain source evidence changed '+str(key))
  if row['source_status']=='CNKI_user_confirmed' or row.get('source_declaration_file'):
   declaration=Path(row['source_declaration_file']);require(declaration.is_file(),'Source declaration missing '+str(key))
   if str(declaration) not in source_hashes:source_hashes[str(declaration)]=sha(declaration)
   require(source_hashes[str(declaration)]==row['source_declaration_sha256'],'Source declaration changed '+str(key))
  budgets.setdefault(y,{})[code]=t
 require({(y,c) for y in required_years for c in CODES}<=seen,'Missing annual provincial budget')
 if required_years==YEARS:require(len(rows)==806 and seen=={(y,c) for y in YEARS for c in CODES},'Formal production requires exactly 806 rows')
 return lock,budgets

def qa_gate(path):
 q=load(path);require(q.get('policy_frozen') is True,'Scientific QA freeze absent')
 require(q.get('qa_rule')=='QA01_DETAILED' and q.get('primary_rule')=='QA01_DETAILED_n10g4' and q.get('grid_min_qualified_native_area_fraction')==.8,'Unsupported QA rule')
 require(q.get('native_min_valid_windows')==10 and q.get('native_min_time_groups')==4 and q.get('primary_ev_band')==1,'Primary temporal rule mismatch')
 require(q.get('sensitivity_rule')=='QA01_DETAILED_n5g3' and q.get('sensitivity_ev_band')==3 and q.get('sensitivity_grid_min_qualified_native_area_fraction')==.8,'Sensitivity rule mismatch')
 require(q.get('composite_start_DOY')==list(range(49,354,16)) and q.get('spatial_support',{}).get('missing_year_imputation') is False,'Temporal support mismatch')
 require(q.get('reference_years')==REFERENCE and q.get('selection_used_EHI_outcomes') is False,'Anchor or rule-selection mismatch')
 require(q.get('FP_missing_weight_model') and q.get('evidence'),'FP/QA decision evidence absent')
 for e in q['evidence']:require(sha(e['path'])==e['sha256'],'QA decision evidence changed')
 return dict(**q,primary_support='n10g4')

def validate_proof(p):
 require(p.get('sha256') and sha(p['path'])==p['sha256'],'Missing/changed receipt '+p['path']);r=load(p['path'])
 for k in p['checks']:require(r.get(k) is True,'Receipt failed '+k+': '+p['path'])
 walk_bindings(r)
def validate_spec(s,grid):
 require(s.get('fingerprint') is not None,'Missing raster '+s['path']);check_fp(s['fingerprint']);require(str(Path(s['path']).resolve())==s['fingerprint']['path'],'Raster fingerprint path differs')
 for p in s['proofs']:validate_proof(p)
 with rasterio.open(s['path']) as src:
  require((src.shape,src.crs,src.transform)==(grid.shape,grid.crs,grid.transform),'Grid mismatch '+s['path']);require(list(src.descriptions)==s['schema'],'Band schema changed '+s['path']);require(src.nodata is not None,'Missing nodata schema')

def contains_fingerprint(obj,expected):
 if isinstance(obj,dict):
  if all(obj.get(k)==v for k,v in expected.items()):return True
  return any(contains_fingerprint(v,expected) for v in obj.values())
 if isinstance(obj,list):return any(contains_fingerprint(v,expected) for v in obj)
 return False

def preflight(config,pilot_year=None):
 cfg=load(config);errors=[];qa=None;budgets=None
 selected=YEARS if pilot_year is None else [pilot_year]
 for name,fun in [('QA',lambda:qa_gate(cfg['qa_lock'])),('grain',lambda:grain_gate(cfg['grain_lock'],selected))]:
  try:
   result=fun()
   if name=='QA':qa=result
   else:_,budgets=result
  except (OSError,ValueError,KeyError,GateError) as e:errors.append(name+': '+str(e))
 try:
  require(cfg['reference_years']==REFERENCE and cfg['formula']==FORMULA,'Formula/reference changed');require(set(map(str,selected))<=set(cfg['years']),'Requested year manifest incomplete')
  if pilot_year is None:
   require(set(cfg['years'])==set(map(str,YEARS)),'26-year manifest incomplete');validate_proof(cfg['collection_proof'])
  require(cfg['fp_method']=='same_province_observed_crop_mean_negative_zero_area_fallback','Unreviewed FP model')
  with rasterio.open(cfg['grid']['path']) as grid:
   require(abs(abs(grid.transform.a*grid.transform.e)-1e6)<1e-6 and grid.transform.b==0 and grid.transform.d==0,'Need 1 km nonrotated equal area grid');require(grid.crs.is_projected and 'Albers' in grid.crs.to_wkt(),'Need known equal-area Albers CRS')
   validate_spec(cfg['grid'],grid);validate_spec(cfg['land2020'],grid)
   for y in map(str,selected):
    entry=cfg['years'][y]
    require(set(entry)=={'EV','EO','HQ','Carbon','SR','WY','crop'},'Incomplete seven-component inputs '+y)
    for k,s in entry.items():
     try:
      validate_spec(s,grid)
      require(s['proofs'],'Missing source validation receipts')
      records=[load(p['path']) for p in s['proofs']]
      require(any(contains_fingerprint(r,s['fingerprint']) for r in records),'Receipts do not bind this output')
      require(all('year' not in r or r['year']==int(y) for r in records),'Receipt year mismatch')
      if k=='EV':require(s['schema']==EV_NAMES,'EV semantic schema')
      if k=='crop':require(s['schema']==CROP_NAMES,'Crop semantic schema')
      if k in PHYSICAL_NAMES:
       n=PHYSICAL_NAMES[k];require(s['schema'][:len(n)]==n and (len(s['schema'])==len(n) or (k=='HQ' and int(y) in REFERENCE+[2025] and len(s['schema'])==3)),'Physical semantic schema '+k)
     except (OSError,ValueError,KeyError,GateError) as e:errors.append(y+' '+k+': '+str(e))
 except (OSError,ValueError,KeyError,GateError) as e:errors.append('manifest: '+str(e))
 return dict(passed=not errors,errors=errors,mode='metadata_only_no_raster_pixels_read',years=selected,grain_rows_required=31*len(selected),formal_outputs_created=False),cfg,qa,budgets

def guard(dest):
 import psutil
 require(shutil.disk_usage(Path(dest).anchor).free>=6*1024**3,'Runtime requires at least 6 GiB free disk')
 require(psutil.virtual_memory().available>=2*1024**3,'Runtime requires 2 GiB available RAM')
def verify_sources(cfg,years=None):
 years=YEARS if years is None else years
 for s in [cfg['grid'],cfg['land2020']]+[s for y in years for s in cfg['years'][str(y)].values()]:check_fp(s['fingerprint'])
def storage_plan(cfg,dest,pilot_year=None):
 selected=YEARS if pilot_year is None else [pilot_year];evidence=[]
 with rasterio.open(cfg['grid']['path']) as ds:cells=ds.width*ds.height
 for y in selected:
  for k,s in cfg['years'][str(y)].items():
   if s.get('fingerprint') and s.get('schema'):
    ratio=s['fingerprint']['bytes']/(cells*len(s['schema'])*4)
    evidence.append(dict(year=y,component=k,bytes=s['fingerprint']['bytes'],bands=len(s['schema']),compressed_to_raw_ratio=ratio))
 ratio=min(1.05,max(.15,max(v['compressed_to_raw_ratio'] for v in evidence)*1.30))
 remaining=[]
 for y in selected:
  if not (Path(dest)/f'FP_{y}.json').exists():remaining.append(dict(name=f'FP_{y}',raw_bytes=cells*11*4))
  if pilot_year is None:
   for s in SUPPORTS:
    if not (Path(dest)/f'EHI_{y}_{s}.json').exists():remaining.append(dict(name=f'EHI_{y}_{s}',raw_bytes=cells*len(EHI_NAMES)*4))
 if pilot_year is None and not (Path(dest)/'complete_series_domains.json').exists():remaining.append(dict(name='complete_series_domains',raw_bytes=cells*3))
 estimated=sum(v['raw_bytes'] for v in remaining)*ratio
 return dict(grid_cells=cells,selected_years=selected,source_compression_evidence=evidence,estimated_ratio_with_30pct_margin=ratio,remaining_outputs=remaining,estimated_additional_bytes=int(estimated),safety_reserve_bytes=6*1024**3,required_start_free_bytes=max(12*1024**3,int(estimated)+6*1024**3),estimate_is_not_guarantee=True)
def start_guard(dest,plan):
 import psutil
 require(shutil.disk_usage(Path(dest).anchor).free>=plan['required_start_free_bytes'],'Start disk estimate plus 6 GiB reserve not met')
 require(psutil.virtual_memory().available>=6*1024**3,'Start requires at least 6 GiB available RAM')
def binding(config,cfg):return dict(config_sha256=sha(config),runner_sha256=sha(__file__),qa_lock_sha256=sha(cfg['qa_lock']),grain_lock_sha256=sha(cfg['grain_lock']))

def allocate_year(y,cfg,budgets,dest,bind):
 stage_event(dest,'FP_start',year=y)
 target=dest/f'FP_{y}.tif';done=target.with_suffix('.json')
 if done.exists():
  r=load(done);require(r['complete'] and r['binding']==bind and r['year']==y,'FP reuse gate');check_fp(r['output']);return target
 require(not target.exists() and not target.with_suffix('.partial.tif').exists(),'Unqualified FP output exists')
 sums=np.zeros((100,8));miss=np.zeros((100,2))
 with rasterio.open(cfg['years'][str(y)]['crop']['path']) as crop,rasterio.open(cfg['grid']['path']) as ids:
  stage_event(dest,'FP_denominator',year=y)
  for w in windows(ids):guard(dest);accumulate_crop(sums,miss,ids.read(1,window=w),crop.read(window=w).astype(float))
  pars=allocation_parameters(sums,miss,budgets[y]);partial=target.with_suffix('.partial.tif')
  stage_event(dest,'FP_write',year=y)
  with rasterio.open(partial,'w',**profile(ids,11)) as out:
   for w in windows(ids):guard(dest);a=allocate_block(ids.read(1,window=w),crop.read(window=w).astype(float),pars);out.write(np.where(np.isfinite(a),a,-9999).astype('float32'),window=w)
   for n,name in enumerate(FP_NAMES,1):out.set_band_description(n,name)
   out.update_tags(year=y,quantity='tonnes per 1km cell',model='same-province observed crop area mean; negative weight zero; missing modeled; explicit area fallback')
  totals=np.zeros((100,5))
  stage_event(dest,'FP_full_read_province_closure',year=y)
  with rasterio.open(partial) as ds:
   for w in windows(ids):
    guard(dest);pid=ids.read(1,window=w);a=ds.read([1,2,3,4,5],window=w).astype(float);a[:,pid==0]=0;require((a>=0).all() and np.isfinite(a).all(),'Invalid persisted FP')
    for b in range(5):totals[:,b]+=np.bincount(pid.ravel(),weights=a[b].ravel(),minlength=100)
  closure=[]
  for c,t in budgets[y].items():
   error=np.abs(totals[c]-t)/max(t,1.);require((error<=1e-6).all(),'Persisted FP budget not closed '+str((y,c)))
   closure.append(dict(province_code=c,budget_tonnes=t,persisted_tonnes=totals[c].tolist(),relative_error=error.tolist()))
 verify_sources(cfg,[y]);atomic_replace(partial,target);dump(done,dict(complete=True,year=y,binding=bind,output=fp(target),province_closure=closure,allocation_diagnostics=pars[3],full_read_passed=True));stage_event(dest,'FP_complete',year=y,output=str(target));return target

def raw_block(src,w,support,inside):
 index=SUPPORTS.index(support);a={}
 for k in KEYS:
  band=(1+2*index if k=='EV' else 1+index if k=='SR' else 2+2*index if k=='FP' else 1)
  v=read(src[k],band,w)
  if k=='FP':v/=100. # tonnes/cell to tonnes/ha, same historical physical unit
  v[~inside]=np.nan;a[k]=v
 return a
def open_year(stack,cfg,y,dest):return {k:stack.enter_context(rasterio.open(dest/f'FP_{y}.tif' if k=='FP' else cfg['years'][str(y)][k]['path'])) for k in KEYS}

def fit_anchors(cfg,qa,dest,bind):
 stage_event(dest,'anchors_start',reference_years=REFERENCE)
 path=dest/'fixed_historical_anchors.json'
 if path.exists():
  r=load(path);require(r['complete'] and r['binding']==bind and r['reference_years']==REFERENCE,'Anchor reuse gate');return r
 anc=empty_anchors()
 with rasterio.open(cfg['grid']['path']) as ids,rasterio.open(cfg['land2020']['path']) as land:
  for y in REFERENCE:
   stage_event(dest,'anchors_reference_year',year=y)
   with ExitStack() as st:
    src=open_year(st,cfg,y,dest)
    for w in windows(ids):guard(dest);inside=(ids.read(1,window=w)>0)&(land.read(1,window=w)>0);update_anchors(anc,raw_block(src,w,qa['primary_support'],inside),inside,y)
 require(all(np.isfinite(v['low']) and v['high']>v['low'] for v in anc.values()),'Degenerate reference anchors')
 r=dict(complete=True,binding=bind,reference_years=REFERENCE,primary_support=qa['primary_support'],domain='province>0 and CLCD2020>0; each component finite independently',anchors=anc,anchor_digest=digest(anc),future_or_2025_used=False)
 verify_sources(cfg);dump(path,r);return r

def assemble_year(y,cfg,dest,anc,bind):
 for support in SUPPORTS:
  stage_event(dest,'EHI_start',year=y,support=support)
  target=dest/f'EHI_{y}_{support}.tif';done=target.with_suffix('.json')
  if done.exists():
   r=load(done);require(r['complete'] and r['binding']==bind and r['anchor_digest']==anc['anchor_digest'],'EHI reuse gate');check_fp(r['output']);continue
  require(not target.exists() and not target.with_suffix('.partial.tif').exists(),'Unqualified EHI output exists')
  partial=target.with_suffix('.partial.tif');expected_hash=hashlib.sha256();counts=0
  with ExitStack() as st:
   src=open_year(st,cfg,y,dest);ids=st.enter_context(rasterio.open(cfg['grid']['path']));land=st.enter_context(rasterio.open(cfg['land2020']['path']))
   with rasterio.open(partial,'w',**profile(ids,len(EHI_NAMES))) as out:
    for w in windows(ids):
     guard(dest);inside=(ids.read(1,window=w)>0)&(land.read(1,window=w)>0);a=score(raw_block(src,w,support,inside),anc['anchors']);a[:,~inside]=np.nan;counts+=int(np.isfinite(a[0]).sum());a=np.where(np.isfinite(a),a,-9999).astype('float32');expected_hash.update(a.tobytes());out.write(a,window=w)
    for b,name in enumerate(EHI_NAMES,1):out.set_band_description(b,name)
    out.update_tags(year=y,support=support,formula=FORMULA,anchor_digest=anc['anchor_digest'],anchor_reference_years=','.join(map(str,REFERENCE)),FP_unit='tonnes/ha',primary_support=anc['primary_support'])
   actual_hash=hashlib.sha256()
   with rasterio.open(partial) as out:
    for w in windows(ids):
     guard(dest);a=out.read(window=w);actual_hash.update(a.tobytes());require(np.isfinite(a).all(),'Nonfinite encoded EHI');good=a[0]!=-9999;require(((a[0,good]>=0)&(a[0,good]<=1)).all(),'EHI range');require(np.allclose(a[0,good],a[1,good]*a[2,good],rtol=2e-7,atol=1e-8),'EHI product identity')
   require(actual_hash.hexdigest()==expected_hash.hexdigest(),'Persisted EHI differs from produced blocks')
  verify_sources(cfg);atomic_replace(partial,target);dump(done,dict(complete=True,year=y,support=support,binding=bind,anchor_digest=anc['anchor_digest'],output=fp(target),full_read_passed=True,valid_cells=counts,decoded_sha256=actual_hash.hexdigest()));stage_event(dest,'EHI_complete',year=y,support=support)

def common_domains_and_means(cfg,dest,bind):
 """Different yearly maps never become an apparent temporal trend in means."""
 stage_event(dest,'common_domains_start');target=dest/'complete_series_domains.tif';done=target.with_suffix('.json')
 require(not target.exists() and not done.exists(),'Domain/statistics already exist: inspect receipt, do not overwrite')
 partial=target.with_suffix('.partial.tif');require(not partial.exists(),'Unqualified domain partial exists')
 counts=np.zeros((100,3),np.int64);h=hashlib.sha256()
 with ExitStack() as st:
  ids=st.enter_context(rasterio.open(cfg['grid']['path']));ds={s:[st.enter_context(rasterio.open(dest/f'EHI_{y}_{s}.tif')) for y in YEARS] for s in SUPPORTS}
  p=profile(ids,3);p.update(dtype='uint8',nodata=255,predictor=1)
  with rasterio.open(partial,'w',**p) as out:
   for w in windows(ids):
    guard(dest);pid=ids.read(1,window=w);m=[]
    for s in SUPPORTS:m.append(np.logical_and.reduce([np.isfinite(read(d,1,w)) for d in ds[s]])&(pid>0))
    m.append(m[0]&m[1]);a=np.stack(m).astype('uint8');h.update(a.tobytes());out.write(a,window=w)
    for b in range(3):counts[:,b]+=np.bincount(pid.ravel(),weights=a[b].ravel(),minlength=100).astype('int64')
   for b,n in enumerate(['n10g4_26year_complete','n5g3_26year_complete','both_rules_26year_complete'],1):out.set_band_description(b,n)
  actual=hashlib.sha256()
  with rasterio.open(partial) as out:
   for w in windows(ids):actual.update(out.read(window=w).tobytes())
  require(actual.hexdigest()==h.hexdigest(),'Persisted common-domain mismatch')
 atomic_replace(partial,target);rows=[];stage_event(dest,'common_domain_means_start')
 with rasterio.open(cfg['grid']['path']) as ids,rasterio.open(target) as domain:
  for y in YEARS:
   for si,s in enumerate(SUPPORTS):
    stats={k:np.zeros((100,2)) for k in ('annual_valid','own_complete_series','both_rules_complete_series')}
    with rasterio.open(dest/f'EHI_{y}_{s}.tif') as ds:
     for w in windows(ids):
      guard(dest);pid=ids.read(1,window=w);v=read(ds,1,w);d=domain.read(window=w);masks=[np.isfinite(v)&(pid>0),d[si]==1,d[2]==1]
      for label,mask in zip(stats,masks):
       require(np.isfinite(v[mask]).all(),'Common-domain contains missing EHI');stats[label][:,0]+=np.bincount(pid[mask],minlength=100);stats[label][:,1]+=np.bincount(pid[mask],weights=v[mask],minlength=100)
    for label,table in stats.items():
     for c in [0]+CODES:
      n,total=table[CODES].sum(axis=0) if c==0 else table[c];rows.append(dict(year=y,support=s,domain=label,province_code=c,valid_cells=int(n),EHI_mean=float(total/n) if n else None))
 path=dest/'annual_EHI_means_by_support_domain.csv'
 with path.open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 dump(done,dict(complete=True,binding=bind,output=fp(target),table=fp(path),domain_cells_by_province={str(c):counts[c].tolist() for c in CODES},decoded_sha256=h.hexdigest(),primary_comparison='n10g4 own_complete_series',missing_imputation=False,inference_pq_performed=False))

def authorization_gate(auth,bind,pilot_year=None):
 require(auth.get('allowed') is True and auth.get('root_reviewed') is True and auth.get('binding')==bind,'Production not separately root authorized')
 require(auth.get('max_workers')==1,'Only one assembly worker authorized')
 if pilot_year is None:require(auth.get('scope')=='annual_FP_anchors_ES4_2000_2025','Full production scope mismatch')
 else:require(auth.get('scope')=='annual_FP_pilot' and auth.get('year')==pilot_year,'Pilot authorization must name this exact year')

def run(config,pilot_year=None):
 report,cfg,qa,budgets=preflight(config,pilot_year);dump(HERE/'last_launch_preflight.json',report);require(report['passed'],'Preflight blocked: '+json.dumps(report['errors'],ensure_ascii=False))
 bind=binding(config,cfg);authpath=cfg['production_authorization'] if pilot_year is None else cfg.get('pilot_authorization',str(HERE/'PILOT_AUTHORIZATION.json'))
 authorization_gate(load(authpath),bind,pilot_year)
 base=Path(cfg['output_root']).resolve();require(base.parent==DATA.resolve() and base.name=='annual_ehi_assembly_v1','Output must be isolated canonical destination')
 dest=base if pilot_year is None else base/'pilots'/str(pilot_year)
 plan=storage_plan(cfg,dest,pilot_year);dump(HERE/'last_storage_plan.json',plan);start_guard(dest,plan);dest.mkdir(parents=True,exist_ok=True);lock=dest/'active.lock'
 with lock.open('x') as f:f.write(str(os.getpid()))
 try:
  stage_event(dest,'started',mode='full' if pilot_year is None else 'FP_pilot',selected_years=YEARS if pilot_year is None else [pilot_year],binding=bind,storage_plan=plan)
  if pilot_year is not None:
   allocate_year(pilot_year,cfg,budgets,dest,bind);require(binding(config,cfg)==bind,'Pilot freeze changed');grain_gate(cfg['grain_lock'],[pilot_year]);qa_gate(cfg['qa_lock'])
   dump(dest/'complete.json',dict(complete=True,year=pilot_year,binding=bind,scope='single-year FP only',anchors_or_EHI_computed=False));stage_event(dest,'pilot_complete',year=pilot_year,anchors_or_EHI_computed=False);return
  for y in YEARS:allocate_year(y,cfg,budgets,dest,bind)
  anchors=fit_anchors(cfg,qa,dest,bind)
  for y in YEARS:assemble_year(y,cfg,dest,anchors,bind)
  domain_done=dest/'complete_series_domains.json'
  if domain_done.exists():
   d=load(domain_done);require(d['complete'] and d['binding']==bind,'Domain reuse gate');check_fp(d['output']);check_fp(d['table'])
  else:common_domains_and_means(cfg,dest,bind)
  require(binding(config,cfg)==bind,'Freeze changed during production');grain_gate(cfg['grain_lock']);qa_gate(cfg['qa_lock']);dump(dest/'complete.json',dict(complete=True,binding=bind,years=YEARS,formula=FORMULA,anchor_digest=anchors['anchor_digest'],scope='annual FP and EHI assembly only; no trend/inference'))
  stage_event(dest,'complete',anchors_or_EHI_computed=True)
 except Exception as e:
  stage_event(dest,'failed_or_resource_stop',error=repr(e));raise
 finally:lock.unlink()

if __name__=='__main__':
 ap=argparse.ArgumentParser();g=ap.add_mutually_exclusive_group(required=True);g.add_argument('--prepare',type=Path);g.add_argument('--preflight',type=Path);g.add_argument('--run',type=Path);g.add_argument('--pilot',type=int,choices=YEARS);ap.add_argument('--config',type=Path);ap.add_argument('--year',type=int,choices=YEARS,help='Scope metadata preflight to a pilot year');ap.add_argument('--report',type=Path);args=ap.parse_args()
 if args.prepare:prepare(args.prepare)
 elif args.preflight:
  report,*_=preflight(args.preflight,args.year)
  if args.report:dump(args.report,report)
  print(json.dumps(report,ensure_ascii=False,indent=2));raise SystemExit(0 if report['passed'] else 2)
 elif args.pilot is not None:
  if not args.config:ap.error('--pilot requires --config')
  run(args.config,args.pilot)
 else:run(args.run)
