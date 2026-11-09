Current step-0 decision: PASS under the user-authorized raster/count/rounded-percentage agreement criterion. The historical strict 1e-6 test and its precision differences are preserved below. Task 2 additionally verified exact original EHI raster agreement in all 26 years.

# Reproduction log

Step 0 strict gate: FAIL. Runtime: 42.92 seconds.
Original normalization, scoring, raw-band selection and mismatch functions were imported without editing. All retained pixels were re-evaluated from the seven physical component rasters for the four required years.

| Metric | User target | Reproduced | Absolute error | Strict pass | Display agreement |
|---|---:|---:|---:|---|---|
| primary_cells | 4072765 | 4072765 | 0 | True | True |
| reference_cells | 9379006 | 9379006 | 0 | True | True |
| EHI_mean_2000 | 0.023471 | 0.0234707764281 | 2.236e-07 | True | True |
| EHI_mean_2020 | 0.030818 | 0.0308182311817 | 2.312e-07 | True | True |
| EHI_mean_2025 | 0.029976 | 0.0299761262704 | 1.263e-07 | True | True |
| change_2020_2025 | -2.73 | -2.73248943572 | 0.002489 | False | True |
| mismatch_count_2020_2025 | 577714 | 577714 | 0 | True | True |
| mismatch_percent_2020_2025 | 14.18 | 14.1848105648 | 0.004811 | False | True |
| mismatch_count_2015_2020 | 66816 | 66816 | 0 | True | True |
| mismatch_percent_2015_2020 | 1.64 | 1.64055623145 | 0.0005562 | False | True |

Full precision is retained here to diagnose the requested 1e-6 tolerance; four significant digits would hide the discrepancy. Percentages are compared in percentage-point units, not silently redefined as fractions.

Invalid pixels: {2000: 0, 2015: 0, 2020: 0, 2025: 0}. Maximum pixel differences from saved original rasters: {2000: 0.0, 2015: 0.0, 2020: 0.0, 2025: 0.0}.

Tasks 1–7: 未完成（未启动）. The user explicitly requires stopping if step 0 fails. No threshold, target, or input has been adjusted.
## User-authorized gate decision (2026-11-09)

Step 0: PASS under the subsequently authorized criterion: exact original raster agreement, exact cell counts, and agreement of rounded percentages. The strict numerical comparison above is retained unchanged. The earlier task-not-started status describes the historical stop, not the current run. Tasks 1–6 are now proceeding in order; no input, threshold or result was corrected to force agreement.


## Completed supplementary analyses

| Task | Status | Runtime (s) | Cells | Invalid or undefined |
|---|---|---:|---:|---:|
| 1 | complete with explicit cartographic boundary | 82.83 | 4072765 | 0 |
| 2 | complete | 415.2 | 4072765 | 0 |
| 3 | complete | 270.2 | 4072765 | 0 |
| 4 | complete | 746.4 | 4072765 | 0 |
| 5 | complete with undefined variance cells reported | 45 | 4072765 | 48375 |
| 6 | complete | 2.562 | 4072765 | 0 |

Task 1 north overlay has a separately recorded runtime in task1_audit.json. All published CSV floating values use four significant digits; exact integer counts and full-precision audit caches are retained. Zero newly invalid primary-domain cells occurred in tasks 2–4 and 6. Undefined Hamed–Rao variances, if present, remain in the denominator and are listed separately, with p=1 used only conservatively in the BH ordering.

The Qinling–Huaihe line is a reproducible cartographic partition with western extension, not an official administrative boundary. The climate-invalid reference cells are reported separately and do not occur in the primary domain. Fixed-climatology task 2 differs from the manuscript’s prior fixed-endpoint-weather comparison and must not replace it without changing the description. Task 4 uses the same full primary domain as the prior WY summary and must not be conflated with the smaller GLDAS-matched subset.

Final category-conservation validation identified one significant MK decline with zero Sen slope because of ties. Classification was corrected to use the MK score direction for this case, matching pyMannKendall; input values, slopes and p values were unchanged. An unsupported superscript font glyph in the trend-map label was also repaired without altering its data.
