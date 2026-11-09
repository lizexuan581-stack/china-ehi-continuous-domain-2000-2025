"""Mask overlays only. Aridity is ratio of 26-year mean P to mean ET0."""
from common import *
t0=time.monotonic()
mask=read(src['common26'],1)>0;prov=read(src['province']);land=read(src['land2020']);ref=(prov>0)&(land>0)
assert mask.sum()==N and ref.sum()==9379006 and not (mask&~ref).any()
np.save(CACHE/'main_mask.npy',mask);np.save(CACHE/'reference_mask.npy',ref)
ps=np.zeros(mask.shape,'float64');es=np.zeros_like(ps);valid=np.ones(mask.shape,bool)
for y in YEARS:
 p,e=read(climatepath(y),1),read(climatepath(y),2);valid&=np.isfinite(p)&np.isfinite(e)&(e>0);ps+=np.nan_to_num(p);es+=np.nan_to_num(e);print('climate',y,flush=True)
ps/=26;es/=26;ps[~valid]=np.nan;es[~valid]=np.nan
np.save(CACHE/'mean_P.npy',ps);np.save(CACHE/'mean_ET0.npy',es)
ai=np.divide(ps,es,out=np.full(ps.shape,np.nan),where=valid)
# UNEP-style four classes: arid includes hyper-arid (<0.05).
zone=np.zeros(mask.shape,'uint8');zone[valid]=np.digitize(ai[valid],[.2,.5,.65])+1
rows=[]
names={int(r['province_code']):r['province_name'] for r in csv.DictReader((ROOT/'metadata/province_names.csv').open(encoding='utf-8-sig'))}
for code in sorted(names):
 if code==0:continue
 r=ref&(prov==code);c=mask&r;rows.append(dict(group='province',region=names[code],reference_cells=int(r.sum()),main_cells=int(c.sum()),coverage_pct=100*c.sum()/r.sum(),share_of_main_pct=100*c.sum()/N,status='complete'))
for z,name in enumerate(['climate_invalid','arid_AI_lt_0.2','semi_arid_0.2_to_0.5','dry_subhumid_0.5_to_0.65','humid_AI_ge_0.65']):
 r=ref&(zone==z);c=mask&(zone==z);rows.append(dict(group='aridity',region=name,reference_cells=int(r.sum()),main_cells=int(c.sum()),coverage_pct=100*c.sum()/r.sum() if r.sum() else 0.,share_of_main_pct=100*c.sum()/N,status='complete'))
csvout('task1_domain_composition',rows);np.save(CACHE/'aridity_zones.npy',zone)
plt=plot_setup();from matplotlib.colors import ListedColormap;from matplotlib.patches import Patch
fig,ax=plt.subplots(figsize=(7.2,6.7));show=np.where(mask,zone,np.where(ref,5,0));ax.imshow(show,cmap=ListedColormap(['white','#d8b365','#f6e8c3','#80cdc1','#35978f','#dedede']),vmin=0,vmax=5,interpolation='nearest');ax.set_axis_off()
ax.legend(handles=[Patch(color=c,label=l) for c,l in zip(['#d8b365','#f6e8c3','#80cdc1','#35978f','#dedede'],['Arid (AI < 0.2)','Semi-arid (0.2–0.5)','Dry subhumid (0.5–0.65)','Humid (AI ≥ 0.65)','Outside continuous observation domain'])],loc='lower left',fontsize=8)
ax.set_title('Observation support by climatic aridity',loc='left');savefig(fig,'task1_domain_composition');plt.close(fig)
record('task1',t0,main_cells=N,reference_cells=int(ref.sum()),climate_invalid_main=int((mask&~valid).sum()),north_south_status='pending sourced boundary',aridity_definition='mean annual P / mean annual ET0, 2000–2025; <0.2, [0.2,0.5), [0.5,0.65), >=0.65')
