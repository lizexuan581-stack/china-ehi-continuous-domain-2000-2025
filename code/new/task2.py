"""Fixed-mean direct-weather experiment; observed annual non-WY inputs retained."""
from common import *
t0=time.monotonic();mask=np.load(CACHE/'main_mask.npy');assert mask.sum()==N
p=np.load(CACHE/'mean_P.npy')[mask].astype('float32');e=np.load(CACHE/'mean_ET0.npy')[mask].astype('float32');st=static(mask)
series=np.lib.format.open_memmap(CACHE/'observed_EHI.npy',mode='w+',dtype='float32',shape=(26,N))
wn=np.lib.format.open_memmap(CACHE/'weather_normalized_EHI.npy',mode='w+',dtype='float32',shape=(26,N))
phis=np.lib.format.open_memmap(CACHE/'PHI.npy',mode='w+',dtype='float64',shape=(26,N))
rows=[];means=[];wm=[];invalid=[];wcoef=[]
for i,y in enumerate(YEARS):
 vals={k:raw(y,k,mask) for k in core.KEYS};scores=np.array([norm(k,vals[k]) for k in core.KEYS]);phi=np.cbrt(scores[0]*scores[1]*scores[2]);phis[i]=phi
 obs=read(DEST/f'EHI_{y}_n10g4.tif')[mask];assert np.array_equal(core.score(scores).astype('float32'),obs)
 series[i]=obs;wcoef.append((phi*scores[3:7]).mean(axis=1))
 newwy=water(read(landpath(y))[mask],p,e,st,1);scores[5]=norm('WY',newwy.astype('float64'));new=core.score(scores).astype('float32');wn[i]=new;bad=int((~np.isfinite(new)).sum());invalid.append(bad)
 if bad:raise RuntimeError('New invalid cells; denominator must not shrink')
 means.append(float(obs.mean(dtype='float64')));wm.append(float(new.mean(dtype='float64')))
 rows.append(dict(record='annual_mean',year=int(y),cells=N,invalid_cells=bad,EHI=means[-1],EHI_wn=wm[-1],WY_wn_mean_mm=float(newwy.mean(dtype='float64'))))
 if y in [2000,2015,2020,2025]:np.savez(CACHE/f'endpoint_{y}.npz',ndvi=vals['EV'],wn_scores=scores)
 print('task2',y,'mean',wm[-1],flush=True)
series.flush();wn.flush();phis.flush();np.save(CACHE/'weight_service_coefficients.npy',wcoef)
for a,b in [(2015,2020),(2020,2025)]:
 na=np.load(CACHE/f'endpoint_{a}.npz');nb=np.load(CACHE/f'endpoint_{b}.npz');mc=mismatch(wn[a-2000].astype(float),wn[b-2000].astype(float),na['ndvi'],nb['ndvi'])
 contrib=core.exact_shapley(na['wn_scores'],nb['wn_scores']);channel=core.negative_channel(contrib)
 for label,nt,et in [('raw',1e-6,1e-6),('screened',.005,.001)]:
  mm=core.classify(nb['ndvi']-na['ndvi'],wn[b-2000].astype(float)-wn[a-2000].astype(float),nt,et)==2
  rows.append(dict(record='mismatch',period=f'{a}–{b}',threshold=label,cells=N,mismatch_cells=mc[label],mismatch_pct=100*mc[label]/N,WY_dominant_cells=int((mm&(channel==1)).sum()),WY_dominant_pct=100*(mm&(channel==1)).sum()/mc[label] if mc[label] else 0.))
 change=wm[b-2000]-wm[a-2000];od=means[b-2000]-means[a-2000]
 rows.append(dict(record='period_change',period=f'{a}–{b}',cells=N,observed_delta=od,wn_delta=change,observed_change_pct=100*od/means[a-2000],wn_change_pct=100*change/wm[a-2000],weather_dependent_difference=od-change,weather_dependent_fraction_of_observed_delta=(od-change)/od))
rows.append(dict(record='regional_Sen',period='2000–2025',cells=N,EHI_slope=sen(means),EHI_wn_slope=sen(wm)))
csvout('task2_weather_normalization',rows)
plt=plot_setup();fig,ax=plt.subplots(figsize=(7.2,4.5));ax.plot(YEARS,means,'o-',ms=3,color='#585858',label='Observed weather');ax.plot(YEARS,wm,'s-',ms=3,color='#007f86',label='Fixed 2000–2025 mean weather');ax.set(xlabel='Year',ylabel='Mean EHI',title='Direct-weather dependence on identical spatial support');ax.legend();fig.tight_layout();savefig(fig,'task2_weather_normalization');plt.close(fig)
record('task2',t0,status='complete',cells=N,invalid_cells=invalid,mean_weather_cast='float32, consistent with original physical-input precision',fixed_anchors=True,causal_attribution=False)
