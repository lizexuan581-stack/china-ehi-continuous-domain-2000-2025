"""Reproduce the approved raster/count/rounded-percentage gate before analysis."""
from common import *
t0=time.monotonic();mask=read(src['common26'])>0;ref=read(absolute('data/masks/reference_domain.tif'))>0
assert mask.sum()==4072765 and ref.sum()==9379006 and not (mask&~ref).any()
means={};ehi={};ndvi={};rows=[]
for y in [2000,2015,2020,2025]:
 vals={k:raw(y,k,mask) for k in core.KEYS};calc=core.score(core.normalize(vals,anc)).astype('float32');saved=read(DEST/f'EHI_{y}_n10g4.tif')[mask]
 assert np.array_equal(calc,saved),f'Original raster mismatch in {y}'
 means[y]=float(calc.mean(dtype='float64'));ehi[y]=calc.astype(float);ndvi[y]=vals['EV'];rows.append(dict(metric=f'EHI_mean_{y}',value=means[y],cells=N,max_raster_difference=0))
for y,target in [(2000,.023471),(2020,.030818),(2025,.029976)]:assert round(means[y],6)==target
change=100*(means[2025]/means[2020]-1);assert round(change,2)==-2.73;rows.append(dict(metric='change_2020_2025_pct',value=change,cells=N))
for a,b,count,pct in [(2015,2020,66816,1.64),(2020,2025,577714,14.18)]:
 n=mismatch(ehi[a],ehi[b],ndvi[a],ndvi[b])['screened'];assert n==count and round(100*n/N,2)==pct
 rows.append(dict(metric=f'mismatch_{a}_{b}',value=100*n/N,mismatch_cells=n,cells=N))
csvout('task0_reproduction',rows);record('task0',t0,status='pass',criterion='Exact original rasters and counts; percentages agree after original display rounding',cells=N,reference_cells=int(ref.sum()))
