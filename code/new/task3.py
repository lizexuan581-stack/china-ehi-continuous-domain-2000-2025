"""FP quantiles fitted on the original component-specific reference support."""
from common import *
t0=time.monotonic();mask=np.load(CACHE/'main_mask.npy');ref=np.load(CACHE/'reference_mask.npy')
pool=[]
for y in [2000,2005,2010,2015,2020]:
 f=read(rawpath(y,'FP'),2).astype('float64')/100;pool.append(f[ref&np.isfinite(f)]);print('FP reference',y,len(pool[-1]),flush=True)
values=np.concatenate(pool);lo,hi=np.quantile(values,[.02,.98],method='linear');refn=len(values);del pool,values
assert hi>lo;bounds={k:dict(v) for k,v in anc.items()};bounds['FP'].update(low=float(lo),high=float(hi))
rows=[dict(record='reference_bounds',reference_years='2000;2005;2010;2015;2020',reference_values=refn,FP_low=float(lo),FP_high=float(hi),unit='t ha-1 of whole grid cell',quantile_method='linear')];means=[];ends={}
for y in YEARS:
 vals={k:raw(y,k,mask) for k in core.KEYS};scores=np.array([norm(k,vals[k],bounds) for k in core.KEYS]);ehi=core.score(scores).astype('float32');bad=int((~np.isfinite(ehi)).sum());assert bad==0;means.append(float(ehi.mean(dtype='float64')))
 rows.append(dict(record='annual_mean',year=int(y),cells=N,invalid_cells=bad,EHI=means[-1]))
 if y in [2020,2025]:ends[y]=(ehi.astype('float64'),vals['EV'])
 print('task3',y,means[-1],flush=True)
mc=mismatch(ends[2020][0],ends[2025][0],ends[2020][1],ends[2025][1])
for key,v in mc.items():rows.append(dict(record='mismatch_2020_2025',threshold=key,cells=N,mismatch_cells=v,mismatch_pct=100*v/N))
rows.append(dict(record='period_change',period='2020–2025',cells=N,delta=means[-1]-means[20],change_pct=100*(means[-1]/means[20]-1)))
csvout('task3_FP_bounds',rows)
base=np.load(CACHE/'observed_EHI.npy',mmap_mode='r').mean(axis=1,dtype='float64');plt=plot_setup();fig,ax=plt.subplots(figsize=(7.2,4.5));ax.plot(YEARS,base,'o-',ms=3,color='#585858',label='Original min–max bounds');ax.plot(YEARS,means,'s-',ms=3,color='#9c755f',label='FP 2nd–98th percentile bounds');ax.set(xlabel='Year',ylabel='Mean EHI',title='Sensitivity to the food-provision reference range');ax.legend();fig.tight_layout();savefig(fig,'task3_FP_bounds');plt.close(fig)
record('task3',t0,status='complete',cells=N,invalid_cells=0,FP_reference_low=lo,FP_reference_high=hi,reference_values=refn,reference_support='Original per-component finite support within province and valid 2020 land mask; not restricted to 26-year main domain')
