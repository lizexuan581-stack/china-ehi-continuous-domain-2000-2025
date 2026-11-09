# Supplementary analysis contract

The original rasters and producer scripts are read-only inputs. New code resides in this directory; outputs reside in the dated reproduction directory. The user approved the revised step-0 gate on 2026-11-09: identical original rasters, identical counts, and matching rounded percentages. The literal 1e-6 comparison remains in the historical log.

## Definitions shared by all tasks

- Main domain: 4,072,765 cells with all 26 annual EHI observations; fixed throughout.
- Reference domain: province ID > 0 and valid 2020 land cover > 0, 9,379,006 cells.
- Equal 1 km² Albers cells: arithmetic means are the original area-weighted means.
- Five reference years: 2000, 2005, 2010, 2015, 2020. Each physical component uses its original finite reference support. A main-domain restriction is not introduced while fitting sensitivity bounds.
- Original EHI algebra, clipping and component treatment are imported from the frozen original modules. EO and HQ retain their original [0,1] treatment.
- FP band 2 is tonnes per whole cell, divided by 100 ha/cell. Its t ha⁻¹ value is not crop yield.
- Raw mismatch thresholds: NDVI 1e-6 and EHI 1e-6. Screened thresholds: NDVI 0.005 and EHI 0.001. Original classification and exact Shapley code are reused.
- Results are computed at full precision. CSV floats have four significant digits, integer counts remain exact. Float32 EHI storage is retained for consistency with the original rasters.

## Figure contracts

All six outputs use Python/matplotlib, 300 dpi PNG, with an editable PDF audit companion. Each is a single scientific panel; legends and colorbars are accessories. There are no biological replicates or inferential error bars in deterministic scenarios.

1. Domain composition: categorical spatial map, showing retained support within the reference mask and aridity categories. AI = mean(P)/mean(ET0), not mean(P/ET0); arid includes hyper-arid cells. The Qinling–Huaihe overlay is a published cartographic line with a western extension and explicit horizontal endpoint closure. The precise north fraction is conditional on this choice.
2. Weather normalization: paired regional annual sequences on identical cells, comparing original weather with fixed 2000–2025 mean P and ET0. The comparison changes only the direct WY inputs; it is not a causal weather attribution.
3. FP bounds: paired regional annual sequences, comparing original bounds with pooled 2nd–98th reference-year percentiles. No interval denotes statistical uncertainty.
4. Z sensitivity: four annual regional sequences, one per Z with separately refitted WY bounds. Original stored WY bands are independently compared against the newly recomputed values.
5. Trend significance: Sen slope raster with significance display. Every cell is tested; display stippling is aggregated into 20 km blocks with a majority of significant observed cells. The color scale saturates at the 98th percentile of absolute slope, with colorbar extensions. No cells are removed from the computation for rendering.
6. Weight perturbation: distributions of changes under 1,000 Dirichlet draws. Box whiskers show 5th–95th percentiles of weight sensitivity, not confidence intervals. Alpha=(5,5,5,5) gives mean weights 0.25 and total concentration 20. Seed=20261109.

## Statistical implementation

Hamed–Rao MK uses Sen detrending, average ranks, tie-corrected MK variance, biased rank autocorrelation, and significant lags among all 1–25 available lags. The accelerated implementation is checked against installed pyMannKendall on 100 reproducibly sampled observed series and tie fixtures. Undefined adjusted variance is explicitly recorded; such cells are not silently classified as nonsignificant. For conservative BH ordering only, undefined p values are represented by 1 while retaining the original domain size.

References: Hamed and Rao (1998), DOI 10.1016/S0022-1694(97)00125-X; pyMannKendall, DOI 10.21105/joss.01556. Aridity category source: IPCC AR6 WGII Figure CCP3.1, https://www.ipcc.ch/report/ar6/wg2/figures/chapter-ccp3/figure-ccp3-001/ . The use of the study's ET0 is explicit.

## Separation from existing manuscript experiments

The new fixed-climatology comparison does not equal the existing fixed-start-weather or fixed-end-weather comparisons. The new Z sensitivity uses the same continuous-observation main domain as the existing full-domain WY summary; GLDAS validation uses a smaller matched subset. Their numbers must be labeled separately, not substituted interchangeably.
