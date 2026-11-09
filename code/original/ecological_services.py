"""Auditable physical-service algebra, without invented biophysical parameters.

Sources: Zhu et al. (2024), doi:10.3390/f15111987, Table 2, equations 6-11;
archived InVEST Annual Water Yield manual, equations 1-3. Callers supply every
regional coefficient. NaN is missing data, never an assumed zero.
"""
import numpy as np


def _arrays(*values):
    return np.broadcast_arrays(*[np.asarray(v, dtype=float) for v in values])


def _nonnegative(*values):
    if any(np.isinf(v).any() or np.any(v[np.isfinite(v)] < 0) for v in values):
        raise ValueError('Inputs must be nonnegative finite numbers or NaN.')


def _fraction(*values):
    _nonnegative(*values)
    if any(np.any(v[np.isfinite(v)] > 1) for v in values):
        raise ValueError('Fraction outside [0,1].')


def carbon_stock_from_fractions(fractions, densities_mg_c_ha, pixel_area_m2):
    """Sum four pools; preserve missing classes, soil depth and unit provenance.

    fractions: (classes, ...), sum=1 on each valid pixel. Densities: (classes, 4),
    order aboveground, belowground, soil, dead. Four finite zeros mean a source
    explicitly assigns zero; NaN means unresolved. Returns stock Mg C and the
    fraction with all four pools supported. Partial sums are NOT total stocks.
    """
    f=np.asarray(fractions,float); d=np.asarray(densities_mg_c_ha,float)
    if f.ndim<2 or d.shape!=(f.shape[0],4):
        raise ValueError('Expected fractions (classes,...) and four carbon pools.')
    _fraction(f); _nonnegative(d)
    area=np.broadcast_to(np.asarray(pixel_area_m2,float),f.shape[1:])
    _nonnegative(area)
    if np.any(area[np.isfinite(area)]<=0):raise ValueError('Pixel area must be positive.')
    complete=np.isfinite(f).all(axis=0)
    if np.any(complete & ~np.isclose(f.sum(axis=0),1,atol=1e-6)):
        raise ValueError('Valid class fractions must sum to one.')
    known=np.isfinite(d).all(axis=1)
    coverage=np.sum(f[known],axis=0)
    density=np.tensordot(np.sum(d[known],axis=1),f[known],axes=(0,0))
    supported=complete & ~np.any(f[~known]>0,axis=0) & np.isfinite(area)
    stock=np.where(supported,density*area/10000,np.nan)
    return stock,np.where(complete,coverage,np.nan)


def soil_retention(r_erosivity, k_erodibility, ls_factor, c_cover, p_practice):
    """Potential and actual sheet/rill erosion, and their difference (t/ha/yr).

    R in MJ mm ha^-1 h^-1 yr^-1; K in t ha h/(ha MJ mm); LS,C,P dimensionless.
    This is RUSLE avoided erosion, not SDR sediment export or trapped sediment.
    """
    r,k,ls,c,p=_arrays(r_erosivity,k_erodibility,ls_factor,c_cover,p_practice)
    _nonnegative(r,k,ls);_fraction(c,p)
    potential=r*k*ls;actual=potential*c*p
    return potential,actual,potential-actual


def water_yield_budyko(precip_mm, et0_mm, kc, pawc, root_depth_mm,
                       restricting_depth_mm, seasonality_z, vegetated):
    """InVEST annual water-yield equations with every calibrated input explicit.

    omega=min(1.25 + Z*AWC/P,5). vegetated is a supplied 0/1 class flag.
    P=0 gives AET=Y=0 under this annual no-storage-change model, not an observed
    drought-water balance. Snow/ice/melt need a separate process model/mask.
    Interpret hydrologically at watershed scale, as the manual specifies.
    """
    p,e,k,pa,rd,sd,z,v=_arrays(precip_mm,et0_mm,kc,pawc,root_depth_mm,
                             restricting_depth_mm,seasonality_z,vegetated)
    _nonnegative(p,e,k,rd,sd,z);_fraction(pa)
    if np.isinf(v).any() or np.any(np.isfinite(v)&~np.isin(v,[0,1])):
        raise ValueError('vegetated flag must be 0, 1 or NaN.')
    if np.any(np.isfinite(z)&(z<=0)):raise ValueError('Seasonality Z must be positive.')
    common=np.isfinite(p)&np.isfinite(e)&np.isfinite(k)&np.isfinite(v)
    veg=common&(v==1)&np.isfinite(pa)&np.isfinite(rd)&np.isfinite(sd)&np.isfinite(z)
    nonveg=common&(v==0)
    aet=np.full(p.shape,np.nan);omega=np.full(p.shape,np.nan)
    pet=e*k
    aet[nonveg]=np.minimum(p[nonveg],pet[nonveg])
    wet=veg&(p>0);dry=veg&(p==0)
    omega[wet]=np.minimum(1.25+z[wet]*np.minimum(rd[wet],sd[wet])*pa[wet]/p[wet],5)
    phi=pet[wet]/p[wet]
    # Stable Fu equation; avoid overflow in phi**omega and cancellation at large phi.
    a=np.zeros(phi.shape);positive=phi>0
    lp=np.log(phi[positive]);w=omega[wet][positive]
    term=np.exp(np.logaddexp(0,w*lp)/w)
    a[positive]=1+phi[positive]-term
    large=positive&(phi>1)
    pl=phi[large];wl=omega[wet][large]
    a[large]=1-pl*np.expm1(np.log1p(np.exp(-wl*np.log(pl)))/wl)
    aet[wet]=p[wet]*np.clip(a,0,np.minimum(1,phi))
    aet[dry]=0
    return p-aet,aet,omega


def allocate_county_food(ndvi, cropland_mask, total_food):
    """Zhu 2024 Eq.11, equal-area cropland cells within ONE county.

    Require all crop NDVI to be finite and nonnegative. Missing or negative
    cropland NDVI stops allocation; no silent clipping or redistribution.
    total_food must be a reported county total for the matching year and unit.
    """
    n=np.asarray(ndvi,float);mask=np.asarray(cropland_mask)
    if n.shape!=mask.shape or mask.dtype!=np.bool_:raise ValueError('Use a same-shape Boolean crop mask.')
    if not np.isfinite(total_food) or total_food<0:raise ValueError('Unknown county total.')
    if not np.isfinite(n[mask]).all() or np.any(n[mask]<0):raise ValueError('Incomplete/negative crop NDVI.')
    denominator=n[mask].sum()
    if denominator<=0:raise ValueError('No positive crop NDVI; allocation is unidentified.')
    out=np.zeros(n.shape);out[mask]=n[mask]/denominator*total_food
    return out


def normalize_fixed_reference(values, low, high, direction):
    """Explicit analysis rule, not a recovered parameter of Zhu 2024.

    Same anchors across years. Returns out-of-reference flag without clipping;
    the EHI aggregator will reject out-of-range scores until policy is resolved.
    """
    if not np.isfinite([low,high]).all() or high<=low:raise ValueError('Invalid fixed anchors.')
    if direction not in ['positive','negative']:raise ValueError('Specify direction explicitly.')
    x=np.asarray(values,float)
    if np.isinf(x).any():raise ValueError('Infinite input.')
    s=(x-low)/(high-low)
    if direction=='negative':s=1-s
    return s,np.isfinite(s)&((s<0)|(s>1))
