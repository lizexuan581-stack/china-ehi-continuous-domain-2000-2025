"""Hamed–Rao MK, all lags, validated against pyMannKendall; global BH FDR."""
from common import *
from numba import njit,prange,set_num_threads
from math import erfc,sqrt
import pymannkendall as mk
set_num_threads(4)

@njit
def one(x):
 n=len(x); dif=np.empty(n*(n-1)//2); k=0; s=0.
 for i in range(n-1):
  for j in range(i+1,n):
   d=x[j]-x[i];dif[k]=d/(j-i);k+=1;s+=np.sign(d)
 slope=np.median(dif); tie=0.
 for i in range(n):
  first=True
  for j in range(i):
   if x[j]==x[i]:first=False
  if first:
   count=0
   for j in range(n):
    if x[j]==x[i]:count+=1
   tie+=count*(count-1)*(2*count+5)
 var=(n*(n-1)*(2*n+5)-tie)/18.
 if var==0:return slope,1.,0.
 detr=x-np.arange(1,n+1)*slope;ranks=np.empty(n)
 for i in range(n):
  less=0;equal=0
  for j in range(n):
   if detr[j]<detr[i]:less+=1
   if detr[j]==detr[i]:equal+=1
  ranks[i]=less+(equal+1)/2.
 ranks-=ranks.mean();den=(ranks*ranks).sum();sni=0.
 if den>0:
  for lag in range(1,n):
   ac=(ranks[:n-lag]*ranks[lag:]).sum()/den
   if abs(ac)>1.959963984540054/sqrt(n):sni+=(n-lag)*(n-lag-1)*(n-lag-2)*ac
 var*=1+2*sni/(n*(n-1)*(n-2))
 if s==0:return slope,1.,var
 if var<=0:return slope,np.nan,var
 z=(s-np.sign(s))/sqrt(var)
 return slope,erfc(abs(z)/sqrt(2)),var

@njit(parallel=True)
def batch(a):
 out=np.empty((a.shape[1],3))
 for j in prange(a.shape[1]):out[j]=one(a[:,j].astype(np.float64))
 return out

t0=time.monotonic();arr=np.load(CACHE/'observed_EHI.npy',mmap_mode='r');assert arr.shape==(26,N)
rng=np.random.default_rng(20261109); checks=[]
for idx in rng.choice(N,100,replace=False):
 x=arr[:,idx].astype(float);a=one(x);b=mk.hamed_rao_modification_test(x)
 assert np.isclose(a[0],b.slope,atol=1e-14,rtol=1e-12)
 assert np.isclose(a[1],b.p,atol=1e-12,rtol=1e-10,equal_nan=True)
 checks.append(abs(a[1]-b.p))
assert one(np.ones(26))[1]==1
for x in [np.repeat(np.arange(13.),2),rng.normal(size=26)]:
 a=one(x);b=mk.hamed_rao_modification_test(x);assert np.isclose(a[1],b.p,atol=1e-12,equal_nan=True)
out=np.lib.format.open_memmap(CACHE/'HR_MK.npy',mode='w+',dtype='float64',shape=(N,3))
for start in range(0,N,100000):
 end=min(start+100000,N);out[start:end]=batch(np.ascontiguousarray(arr[:,start:end]));print('task5 cells',end,flush=True)
out.flush();p=np.asarray(out[:,1]);valid=np.isfinite(p);bad=int((~valid).sum())
# Undefined adjusted variances retain their original cells; p=1 only for BH's denominator.
sortedp=np.sort(np.where(valid,p,1.));passed=np.where(sortedp<=.05*np.arange(1,N+1)/N)[0];cut=float(sortedp[passed[-1]]) if len(passed) else -1.
sig=valid&(p<=cut);up=sig&(out[:,0]>0);down=sig&(out[:,0]<0);non=valid&~sig
# With many tied zeros, a significant MK score can coexist with a zero median slope.
# Classify these cells by MK score, as pyMannKendall does, without changing p or slope.
zero_sen_sig=np.flatnonzero(sig&(out[:,0]==0))
for j in zero_sen_sig:
 x=arr[:,j].astype(float);score=sum(np.sign(x[i+1:]-x[i]).sum() for i in range(25));up[j]=score>0;down[j]=score<0
assert int(up.sum()+down.sum()+non.sum()+(~valid).sum())==N
rows=[]
for name,v in [('significant_increase',up),('significant_decrease',down),('not_significant',non),('undefined_adjusted_variance',~valid)]:
 rows.append(dict(record='pixel_category',category=name,cells=int(v.sum()),denominator=N,share_pct=100*v.sum()/N,FDR_q=.05,BH_p_cutoff=cut))
means=arr.mean(axis=1,dtype='float64');regional=mk.hamed_rao_modification_test(means)
rows.append(dict(record='regional_test',cells=N,period='2000–2025',p=regional.p,Sen_slope=regional.slope,MK_S=regional.s,adjusted_variance=regional.var_s))
csvout('task5_trend_significance',rows);np.save(CACHE/'HR_significant.npy',sig)
mask=np.load(CACHE/'main_mask.npy');grid=np.full(mask.shape,np.nan,'float32');grid[mask]=out[:,0]*1e3
plt=plot_setup();fig,ax=plt.subplots(figsize=(7.2,6.7));lim=float(np.nanpercentile(abs(grid),98));im=ax.imshow(grid,cmap='RdBu',vmin=-lim,vmax=lim,interpolation='nearest');sign_grid=np.zeros(mask.shape,'float32');sign_grid[mask]=sig
# Display-only 20 km stipples mark blocks with >=50% significant observed cells.
rr=[];cc=[]
for r in range(0,mask.shape[0],20):
 for c in range(0,mask.shape[1],20):
  count=mask[r:r+20,c:c+20].sum()
  if count and sign_grid[r:r+20,c:c+20].sum()/count>=.5:
   br,bc=np.where(mask[r:r+20,c:c+20]);j=np.argmin((br-9.5)**2+(bc-9.5)**2);rr.append(r+br[j]);cc.append(c+bc[j])
ax.scatter(cc,rr,s=.3,c='black',linewidths=0,alpha=.65);ax.set_axis_off();ax.set_title('Sen slope and Hamed–Rao MK significance',loc='left');fig.colorbar(im,ax=ax,shrink=.65,label=r'Sen slope ($10^{-3}$ EHI year$^{-1}$)',extend='both');fig.text(.12,.06,'Dots: ≥50% of observed cells significant within a 20 km display block; BH q = 0.05.',fontsize=7);savefig(fig,'task5_trend_significance');plt.close(fig)
record('task5',t0,status='complete' if not bad else 'complete with undefined variance cells reported',cells=N,invalid_test_cells=bad,validation_max_p_difference=float(np.nanmax(checks)),validation_sample=100,pymannkendall_version=mk.__version__,lags='all 1–25',regional_p=regional.p,BH_cutoff=cut,significant_zero_Sen_cells=len(zero_sen_sig),classification='MK score direction; zero-Sen significant ties explicitly retained',significance_display='20 km block majority; all pixel tests retained')
