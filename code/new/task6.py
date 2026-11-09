"""Dirichlet ES4 weights; exact spatial aggregation by linearity."""
from common import *
t0=time.monotonic();seed=20261109;coef=np.load(CACHE/'weight_service_coefficients.npy');assert coef.shape==(26,4)
observed=np.load(CACHE/'observed_EHI.npy',mmap_mode='r').mean(axis=1,dtype='float64');base=coef@np.full(4,.25);err=float(np.max(abs(base-observed)));assert err<1e-8
weights=np.random.default_rng(seed).dirichlet(np.full(4,5.),size=1000);means=weights@coef.T
assert np.isfinite(means).all();rows=[];distributions=[]
for a,b in [(2000,2025),(2020,2025)]:
 change=100*(means[:,b-2000]-means[:,a-2000])/means[:,a-2000];distributions.append(change)
 for i in range(1000):rows.append(dict(record='draw',period=f'{a}–{b}',draw=i+1,cells=N,invalid_cells=0,weight_C=weights[i,0],weight_SR=weights[i,1],weight_WY=weights[i,2],weight_FP=weights[i,3],change_pct=change[i]))
 lo,hi=np.quantile(change,[.05,.95]);rows.append(dict(record='summary',period=f'{a}–{b}',cells=N,invalid_cells=0,draws=1000,seed=seed,total_concentration=20,negative_draws=int((change<0).sum()),negative_fraction=float((change<0).mean()),change_p05_pct=lo,change_p95_pct=hi))
csvout('task6_weight_perturbation',rows)
plt=plot_setup();fig,ax=plt.subplots(figsize=(7.2,4.5));ax.boxplot(distributions,vert=False,whis=(5,95),showfliers=False,tick_labels=['2000–2025','2020–2025']);ax.axvline(0,color='#777777',lw=.8,ls='--');ax.set(xlabel='Change in regional mean EHI (%)',title='Sensitivity to ES4 weights: 1,000 Dirichlet draws');fig.text(.15,.02,'Whiskers: 5th–95th percentiles across weight draws; not confidence intervals.',fontsize=8);fig.tight_layout(rect=(0,.04,1,1));savefig(fig,'task6_weight_perturbation');plt.close(fig)
record('task6',t0,status='complete',cells=N,invalid_cells=0,seed=seed,draws=1000,alpha=[5,5,5,5],concentration_convention='sum(alpha) = 20',equal_weight_mean_max_difference=err,aggregation='exact by linearity of the weighted ES4 mean; double precision before raster storage')
