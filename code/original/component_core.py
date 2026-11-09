"""Standalone seven-component accounting, no legacy imports or implicit sources."""
import numpy as np
KEYS=('EV','EO','HQ','Carbon','SR','WY','FP')
NORMAL={'EV','Carbon','SR','WY','FP'}
def normalize(raw,anchors):
 out=[]
 for k in KEYS:
  v=np.asarray(raw[k],dtype='float64')
  if not np.isfinite(v).all():raise ValueError('Missing supported input '+k)
  if k in NORMAL:
   lo,hi=anchors[k]['low'],anchors[k]['high']
   if not np.isfinite([lo,hi]).all() or hi<=lo:raise ValueError('Invalid anchor')
   v=(v-lo)/(hi-lo)
  out.append(np.clip(v,0,1))
 return np.array(out)
def score(x):
 x=np.asarray(x,dtype='float64')
 if x.shape[0]!=7 or not np.isfinite(x).all() or np.any((x<0)|(x>1)):raise ValueError('Seven finite normalized components required')
 return np.cbrt(x[0]*x[1]*x[2])*np.sum(x[3:7],axis=0)/4
def exact_shapley(a,b):
 # Additivity: four monomials each with four non-dummy players.
 score(a);score(b)
 aa=np.array(a,dtype='float64',copy=True);bb=np.array(b,dtype='float64',copy=True)
 aa[:3]=np.cbrt(aa[:3]);bb[:3]=np.cbrt(bb[:3]);out=np.zeros_like(aa)
 for service in range(3,7):
  members=(0,1,2,service)
  for i in members:
   j,k,l=[x for x in members if x!=i]
   average=(aa[j]*aa[k]*aa[l]+bb[j]*bb[k]*bb[l])/4
   average+=(bb[j]*aa[k]*aa[l]+aa[j]*bb[k]*aa[l]+aa[j]*aa[k]*bb[l]+bb[j]*bb[k]*aa[l]+bb[j]*aa[k]*bb[l]+aa[j]*bb[k]*bb[l])/12
   out[i]+=(bb[i]-aa[i])*average/4
 return out
def classify(ndvi_delta,ehi_delta,ndvi_threshold,ehi_threshold):
 n=np.asarray(ndvi_delta);e=np.asarray(ehi_delta)
 if n.shape!=e.shape or not np.isfinite(n).all() or not np.isfinite(e).all():raise ValueError('Missing response cannot shrink support')
 out=np.zeros(n.shape,'uint8') # either axis near zero, including exact threshold
 up=n>ndvi_threshold;down=n< -ndvi_threshold;rise=e>ehi_threshold;fall=e< -ehi_threshold
 out[up&rise]=1;out[up&fall]=2;out[down&rise]=3;out[down&fall]=4
 return out
def negative_channel(c,tol=1e-12):
 # Disjoint signed groups; water-yield, structural, other-services; mixed=0.
 groups=np.array([c[5],np.sum(c[:3],axis=0),c[3]+c[4]+c[6]])
 neg=np.maximum(-groups,0);total=neg.sum(axis=0);largest=neg.max(axis=0)
 tie=np.sum(abs(neg-largest)<=tol,axis=0)>1
 eligible=(total>0)&(largest>=.5*total)&~tie
 return np.where(eligible,neg.argmax(axis=0)+1,0).astype('uint8')
