# Simulation study

Original simulation, matching and inference code in scientific execution order. No saved data or results are distributed.

The scripts follow the calculation order:

1. `01_network.R`: stochastic-block network, n=300, seed 2026; writes ignored `data/net.RData` and network PDF.
2. `02_covariance.R`: original three dependence covariance matrices; appends them to `data/net.RData`.
3. `03_generate_data.R`: original generator, seeds and treatment-count constraints; defaults to null, 500 datasets per dependence setting. Optional `power` includes original effects 0, 0.1, 0.2, 0.3, 0.4. Stored LM diagnostics are not the null-table inference.
4. `04_matching.R`: original covariate-only (`Full_1`), single-penalty (`Full_2`) and dual-penalty (`P_2`) matching; writes `data/matched_multiprox_bz_0_v*.RData`. Gurobi and its license are required. The original parallel batch and retry settings remain; review `N_CORES` before launching a costly run.
5. `05_prepare_inputs.R`: converts original-format network and matched objects into `inputs/null_dep*.rds`.
6. `06_run_null_inference.R`: the accepted inference runner with only local path/CLI and source-identity changes. Retains `Y~Z`, `Y~Z+factor(subclass)`, naive RI, normal/empirical Oracle p-values, seeds, 500 draws, covariance memoization and the original eta/rho grid.
7. `07_summarize.R`: validates 150 task outputs, summarizes `mean(p < .05)` (not mean p-values), and draws the original sensitivity figure. Writes only ignored RDS/PDF/PNG outputs under `results/`. `LM_Crude` is repeated across designs in the summary; use one copy for the unmatched table row.
8. `08_export_null_grid.R`: deterministic saved-grid export to `data/ri_grid_outputs/`; no inference recomputation.
9. `09_select_optimal_params.R`: original alpha=.05 selection. Independent fixes (eta,rho)=(0,0); Moderate/Strong minimize distance to alpha, retaining the original `slice_min(..., with_ties=FALSE)` ordering.
10. `10_run_ri_power.R`: original power loop, five effects and three dependence settings/designs; writes `data/power_outputs/ri_power_results.csv`.
11. `11_plot_power.R`: original power figures under `data/power_outputs/`.
12. `12_matching_evaluation.R`: original balance, Mahalanobis similarity and network-distance metrics/plots for null matched files; caches and primary figures under `data/diagnosis/`.
13. `13_replot_matching_labels.R`: original manuscript presentation from step 12 caches; writes `cov_mahal.pdf`, `net_avg.pdf`, `matching_eval_all_combined.pdf` under `data/manuscript_figures/`. These names are used by the manuscript. Optional first argument changes only the local figure directory, never manuscript files.

## Commands

Run from this folder. With saved accepted null inputs under `inputs/`, one-dataset inference smoke:

```text
Rscript 06_run_null_inference.R 1 smoke
```

Task IDs 1–50, 51–100 and 101–150 select Independent, Moderate and Strong respectively, ten datasets per task. To regenerate null inference deliberately, execute each ID once (without `smoke`), then summarize. For example, task 1:

```text
Rscript 06_run_null_inference.R 1
# After all 150 tasks exist:
Rscript 07_summarize.R
```

To generate new null inputs, run `Rscript 01_network.R`, then 02, 03, 04 and 05 in order. Both 03 and 04 default to null; both accept `power` to expose all original effects. Gurobi generation/matching is a substantial computation.

Complete power route, after all accepted null tasks are under `tasks/` (or after a deliberate complete null run):

```text
Rscript 08_export_null_grid.R
Rscript 09_select_optimal_params.R
# Network/covariance from steps 01 and 02 must exist in data/net.RData.
Rscript 03_generate_data.R power
Rscript 04_matching.R power
Rscript 10_run_ri_power.R
Rscript 11_plot_power.R
```

Step 03 with `power` includes and overwrites local generated null files: preserve existing runs first. Step 04 skips existing matching outputs, so keep inputs/matches consistent. Step 10 reads original-format matched RData files directly (05 is null-only), stops if its output CSV exists, and skips missing combinations. Inspect its logs/counts before treating plots as complete.

Matching evaluation and manuscript presentation, after null matching:

```text
Rscript 12_matching_evaluation.R
Rscript 13_replot_matching_labels.R
```

Step 12 reuses all five cached tables if present; use caches from the same matching run or move old caches aside before intentional recomputation.

The grid bridge maps `dep_index`, `sim_id`, `eta`, `rho` unchanged, `design` to `subclass`, `p_value` to `p_val`, and sets `bz=0`, `method=RI_Decay_dmin`. CSV names are `ri_grid_results_bz0_v1.csv`, `ri_grid_results_bz0_v2.csv`, and legacy v3 `ri_grid_results.csv`. Only selector-required columns are exported; unavailable raw statistics are not invented. Accepted grid tasks contain P_2 only. Consequently the selector's auxiliary `ri_table2_ri_decay.csv` contains only available Moderate/Strong P_2 rows, not a complete three-design table. Parameter selection for the original power runner is fully supported. Generated CSVs remain ignored.

## Saved-input schema and local placement

Each `inputs/null_depD.rds` is a list with `dep_index` (D in 1:3), `datasets` (500 data frames, each 300 rows, row names `1:300`, no missing values), `covariance` (300-by-300 numeric matrix) and `network_distance` (aligned 300-by-300 matrix, infinite distances permitted). Dataset columns include `Y`, binary `Z`, `X1:X3`, and `subclass_full1`, `subclass_full2`, `subclass_p_2`; each matched set has treated and control units. Y has no ties. Row/matrix order must agree. Keep saved inputs, intermediate data, `tasks/`, `smoke/` and `results/` ignored.

Original author-held inputs are under `_local/oracle-table-figure/inputs/` relative to the repository root. Copy to `simulation/inputs/` only for an intentional local run. Accepted `task_001.rds` through `task_150.rds` may instead be copied from `_local/oracle-table-figure/tasks/` to `simulation/tasks/` for summary-only execution. The three fields used from each task object are `identity`, `table`, `grid`; identity records task/dependence/replications/MC count/seed and hashes. The runner checks its own checkpoint identity exactly. Older accepted task checkpoints remain valid for summary-only use but cannot be resumed by the renamed runner because source-identity paths differ; do not overwrite historical evidence or mix old and new runs.

## Short provenance / read order

- 01–03 map respectively to `code/All in use/simulation/net.R`, `variance.R`, `dga.R`, preserved in `_local/initial-wide-staging-2026-09-29/code/simulation/`. Changes: relative paths, ignored plot destination, explicit covariance load, null-default/power selector, removal of private commented paths. Generator/covariance functions are unchanged.
- 04 maps to `code/All in use/matching/run_matching_batch.R`; `source/utility.R`, `fmnp_gurobi.R`, `run_matching_task.R` are byte-preserved from the matching directory in that same staged lineage. Changes to the batch entrypoint: relative paths and null selector only. Exact archival mapping is recorded in private `_local/PREP_MANIFEST.csv`.
- 05–07 map to accepted `_local/oracle-table-figure/prepare_inputs.R`, `run_null_inference.R`, `collect_compare_plot.R`. 07 retains task validation, rejection-rate calculations and plotting; removes historical CSV comparison/report generation and substitutes ignored RDS output for CSV.
- `source/statistics.R`, `covariance.R`, `utility_oracle.R` are byte-preserved from the accepted reproduction's `source/`. Read them with 06: observed statistic, covariance bound/decay, Oracle draws. The author's accepted evidence remains unchanged under `_local/oracle-table-figure/results/`; no historical-number adjudication is required to use this code.

- 08 is a packaging-only schema bridge. 09–11 map to original `analysis/select_optimal_params.R`, `run_ri_power.R`, `plot_power.R`; 12 to `evaluation/run_matching_eval.R`; 13 to `figures/replot_matching_labels.R`. Changes are limited to paths/source locations and figure destination setup.
- Added `source/utility_ri.R` and `matching_metrics.R` are byte-identical originals. Existing `statistics.R` and `covariance.R` were compared to the original analysis copies and are byte-identical: no duplication is needed.

Dependencies and the verification boundary are in the root README.
