# Ecosystem health on continuously observed land in China, 2000–2025

This repository accompanies a 1 km Albers equal-area assessment of ecosystem health on **4,072,765 continuously observed cells**, representing **43.42% of a 9,379,006-cell analytical reference mask**. Coverage is geographically selective. Results do not describe unobserved southern land or all land in China.

The released data are the **original aligned analytical rasters used in the calculations**, not a claim to include every upstream satellite scene. EHI = PHI × ES4, PHI = (EV × EO × HQ)^(1/3), and ES4 is the equal-weight mean of normalized carbon storage, soil retention, water yield and food provision. Fixed component reference bounds pool 2000, 2005, 2010, 2015 and 2020; 2025 is excluded from fitting. Values outside the bounds are clipped to [0,1]. EO and HQ retain their native scores.

## Data and file layout

Download all `annual_2000.zip` through `annual_2025.zip` and `static_masks_boundaries.zip` from the versioned GitHub release. Extract them **into the repository root**, preserving their `data/` directories. Each release asset is below 2 GiB. `metadata/manifest.json` records every analytical input's path, byte size, SHA-256 and role; the release checksum file identifies the archives.

Each annual archive contains EHI, the seven original component rasters, P/ET0 and land cover. Original multiband GeoTIFFs are retained. EHI is band 1; component inputs use band 1 except FP, which uses band 2 divided by 100. P and ET0 are climate bands 1 and 2. WY bands 1–4 correspond to Z = 1, 5, 10, 20. FP is **tonnes per hectare of whole grid-cell area, not crop yield**.

Static inputs include the original soil and terrain files, class-based root-depth and Kc parameters, the full complete-series domain raster, standalone main/reference masks, and county reporting boundaries. The county source is geoBoundaries CHN ADM3 (2017 representation), licensed ODbL 1.0; it is not an official Chinese administrative boundary product.

The GADM provincial vector and its administrative-ID raster are **not redistributed**. Their official download URL, source checksum and exact reconstruction are provided in `code/new/prepare_province.py`. This respects the original source restriction while allowing permitted academic reproduction. The Qinling–Huaihe overlay uses a published cartographic line with a western extension and documented endpoint closure; its north/south percentage is conditional on that operational partition.

## Reproduction

Install Python 3.10+ and the packages in `requirements.txt`. Original producer scripts are retained without changes under `code/original/`; some retain historical local paths and must not be launched as an unattended legacy production pipeline. Portable analysis scripts under `code/new/` adapt file discovery and figure export while retaining the scientific functions and hash-bound WY helper.

Run these commands from the repository root after extracting the data:

```text
python code/new/task0.py
python code/new/prepare_province.py
python code/new/task1.py
python code/new/task1_north.py
python code/new/task2.py
python code/new/task3.py
python code/new/task4.py
python code/new/task5.py
python code/new/task6.py
```

Stop if task 0 fails. The approved gate requires exact original-raster agreement, exact cell counts and agreement of percentages after their original display rounding. The earlier literal 1e-6 comparison and precision differences remain documented in `reports/reproduction_log.md`.

The tasks cover domain composition, fixed-climatology WY inputs, FP 2nd–98th percentile bounds, Z-specific WY bounds, Hamed–Rao MK with BH FDR, and 1,000 Dirichlet service-weight draws. Weight alpha=(5,5,5,5), total concentration=20, seed=20261109. Equal-area cell means are used throughout. All statistical cells are retained in denominators; 48,375 cells have undefined Hamed–Rao corrected variances and are explicitly reported separately. One significant MK decline has zero Sen slope because of ties and is classified using the MK score direction.

Outputs are CSV and 300 dpi PNG. Published results and audit summaries are in `results/` and `reports/`; fresh runs write to `results_run/`. Full computation caches require additional disk space. Zhu et al.'s county EHI values were not available, so the optional external county correlation was not performed.

## Interpretation and provenance

Fixed-weather comparisons describe propagation through the specified WY calculation; observed vegetation and food maps can retain weather responses. They do not identify causal climate effects. Spatial block composition ranges and weight-perturbation percentiles are not sampling confidence intervals. The new fixed-2000–2025-climatology comparison is distinct from the manuscript's fixed-endpoint-weather experiments.

Primary source links:

- GADM 4.1: https://gadm.org/license.html and https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_CHN_1.json.zip
- County boundary: https://github.com/wmgeolab/geoBoundaries/tree/9469f09/releaseData/gbOpen/CHN/ADM3 ; ODbL 1.0: https://opendatacommons.org/licenses/odbl/1-0/
- Qinling–Huaihe cartographic source: https://github.com/Wsf921109/stata-prov-map ; source license accompanies `code/new/sources/`.
- Hamed–Rao: https://doi.org/10.1016/S0022-1694(97)00125-X ; pyMannKendall: https://doi.org/10.21105/joss.01556
- Aridity categories: https://www.ipcc.ch/report/ar6/wg2/figures/chapter-ccp3/figure-ccp3-001/

The complete historical data and accompanying code/results are archived on Zenodo: **[10.5281/zenodo.24073205](https://doi.org/10.5281/zenodo.24073205)**. All 29 archival files were checked against local checksums before publication. The GitHub release contains 27 data archives and their SHA-256 list; Zenodo additionally contains the code/results ZIP. Third-party licence exclusions described above apply to both copies. The earlier Qinghai-only record is a different study.
