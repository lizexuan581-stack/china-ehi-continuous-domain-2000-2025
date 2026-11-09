"""Recompute WY with original helper for every Z; fit Z-specific WY anchors."""
from common import *
t0=time.monotonic();mask=np.load(CACHE/'main_mask.npy');ref=np.load(CACHE/'reference_mask.npy');zs=[1,5,10,20]
bounds={z:{'low':np.inf,'high':-np.inf,'count':0} for z in zs};stref=static(ref)
for y in [2000,2005,2010,2015,2020]:
 lc=read(landpath(y))[ref];p=read(climatepath(y),1)[ref];e=read(climatepath(y),2)[ref]
 for z in zs:
  wy=water(lc,p,e,stref,z);good=wy[np.isfinite(wy)];b=bounds[z];b['low']=min(b['low'],float(good.min()));b['high']=max(b['high'],float(good.max()));b['count']+=len(good)
 print('task4 fit',y,bounds,flush=True)
del stref;st=static(mask);means={z:[] for z in zs};ends={};maxerrors={z:0. for z in zs};rows=[]
for z,b in bounds.items():rows.append(dict(record='WY_reference_bounds',Z=z,reference_values=b['count'],low_mm=b['low'],high_mm=b['high']))
for y in YEARS:
 vals={k:raw(y,k,mask) for k in core.KEYS};sc=np.array([norm(k,vals[k]) for k in core.KEYS]);lc=read(landpath(y))[mask];p=read(climatepath(y),1)[mask];e=read(climatepath(y),2)[mask]
 for zi,z in enumerate(zs):
  wy=water(lc,p,e,st,z);old=read(rawpath(y,'WY'),zi+1)[mask];err=float(np.max(np.abs(wy.astype(float)-old)));maxerrors[z]=max(maxerrors[z],err)
  # Do not silently accept physical discrepancies from the frozen producer.
  if err>1e-6:raise RuntimeError(f'WY reproduction mismatch: year={y}, Z={z}, error={err}')
  sc[5]=np.clip((wy.astype(float)-bounds[z]['low'])/(bounds[z]['high']-bounds[z]['low']),0,1);ehi=core.score(sc).astype('float32');bad=int((~np.isfinite(ehi)).sum());assert bad==0
  means[z].append(float(ehi.mean(dtype='float64')));rows.append(dict(record='annual_mean',year=int(y),Z=z,cells=N,invalid_cells=bad,EHI=means[z][-1],WY_mm=float(wy.mean(dtype='float64'))))
  if y in [2020,2025]:ends[z,y]=ehi.astype('float64')
 if y in [2020,2025]:ends['ndvi',y]=vals['EV']
 print('task4',y,{z:means[z][-1] for z in zs},flush=True)
for z in zs:
 mc=mismatch(ends[z,2020],ends[z,2025],ends['ndvi',2020],ends['ndvi',2025])
 rows.append(dict(record='period_change',period='2020–2025',Z=z,cells=N,delta=means[z][-1]-means[z][20],change_pct=100*(means[z][-1]/means[z][20]-1)))
 for key,v in mc.items():rows.append(dict(record='mismatch_2020_2025',Z=z,threshold=key,cells=N,mismatch_cells=v,mismatch_pct=100*v/N))
csvout('task4_Z_sensitivity',rows)
plt=plot_setup();fig,ax=plt.subplots(figsize=(7.2,4.5))
for z,c in zip(zs,['#222222','#007f86','#c58c2b','#8663a1']):ax.plot(YEARS,means[z],label=f'Z = {z}',color=c,lw=1.4)
ax.set(xlabel='Year',ylabel='Mean EHI',title='Water-yield parameter sensitivity with refitted WY bounds');ax.legend(ncol=2);fig.tight_layout();savefig(fig,'task4_Z_sensitivity');plt.close(fig)
record('task4',t0,status='complete',cells=N,invalid_cells=0,max_WY_reproduction_error=maxerrors,bounds=bounds,other_component_anchors_unchanged=True)
