# Simulation study

Null-only simulation code in scientific execution order. Accepted OSCAR inference uses the original saved matched datasets, not the superseded generic package runner. No saved data or results are distributed.

The scripts follow the calculation order:

1. `01_network.R`: stochastic-block network, n=300, seed 2026; writes ignored `data/net.RData` and network PDF.
2. `02_covariance.R`: original three dependence covariance matrices; appends them to `data/net.RData`.
3. `03_generate_data.R`: original generator, seeds and treatment-count constraints; null selector; 500 datasets per dependence setting in `data/bz0_v*.RData`. Its stored covariate-adjusted LM diagnostics are historical generator diagnostics, not the null-table inference.
4. `04_matching.R`: original covariate-only (`Full_1`), single-penalty (`Full_2`) and dual-penalty (`P_2`) matching; writes `data/matched_multiprox_bz_0_v*.RData`. Gurobi and its license are required. The original parallel batch and retry settings remain; review `N_CORES` before launching a costly run.
5. `05_prepare_inputs.R`: converts original-format network and matched objects into `inputs/null_dep*.rds`.
6. `06_run_null_inference.R`: the accepted inference runner with only local path/CLI and source-identity changes. Retains `Y~Z`, `Y~Z+factor(subclass)`, naive RI, normal/empirical Oracle p-values, seeds, 500 draws, covariance memoization and the original eta/rho grid.
7. `07_summarize.R`: validates 150 task outputs, summarizes `mean(p < .05)` (not mean p-values), and draws the original sensitivity figure. Writes only ignored RDS/PDF/PNG outputs under `results/`. `LM_Crude` is repeated across designs in the summary; use one copy for the unmatched table row.

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

To generate new inputs instead, run `Rscript 01_network.R`, then 02, 03, 04 and 05 in order. This route is included for inspectability and requires substantial computation; it was not run in full during packaging and is not asserted to regenerate historical solver/Monte Carlo bits. Original generator/batch selectors were changed from nonzero effects to beta_z=0 for the accepted scope. No nonzero-effect run is requested or validated here.

## Saved-input schema and local placement

Each `inputs/null_depD.rds` is a list with `dep_index` (D in 1:3), `datasets` (500 data frames, each 300 rows, row names `1:300`, no missing values), `covariance` (300-by-300 numeric matrix) and `network_distance` (aligned 300-by-300 matrix, infinite distances permitted). Dataset columns include `Y`, binary `Z`, `X1:X3`, and `subclass_full1`, `subclass_full2`, `subclass_p_2`; each matched set has treated and control units. Y has no ties. Row/matrix order must agree. Keep saved inputs, intermediate data, `tasks/`, `smoke/` and `results/` ignored.

Original author-held inputs are under `_local/oracle-table-figure/inputs/` relative to the repository root. Copy to `simulation/inputs/` only for an intentional local run. Accepted `task_001.rds` through `task_150.rds` may instead be copied from `_local/oracle-table-figure/tasks/` to `simulation/tasks/` for summary-only execution. The three fields used from each task object are `identity`, `table`, `grid`; identity records task/dependence/replications/MC count/seed and hashes. The runner checks its own checkpoint identity exactly. Older accepted task checkpoints remain valid for summary-only use but cannot be resumed by the renamed runner because source-identity paths differ; do not overwrite historical evidence or mix old and new runs.

## Short provenance / read order

- 01–03 map respectively to `code/All in use/simulation/net.R`, `variance.R`, `dga.R` in the original MatchingNet source, preserved in `_local/initial-wide-staging-2026-09-29/code/simulation/`. Changes: relative local paths, ignored plot destination, explicit covariance load, null-only selector, removal of private commented paths. Generator/covariance functions are unchanged.
- 04 maps to `code/All in use/matching/run_matching_batch.R`; `source/utility.R`, `fmnp_gurobi.R`, `run_matching_task.R` are byte-preserved from the matching directory in that same staged lineage. Changes to the batch entrypoint: relative paths and null selector only. Exact archival mapping is recorded in private `_local/PREP_MANIFEST.csv`.
- 05–07 map to accepted `_local/oracle-table-figure/prepare_inputs.R`, `run_null_inference.R`, `collect_compare_plot.R`. 07 retains task validation, rejection-rate calculations and plotting; removes historical CSV comparison/report generation and substitutes ignored RDS output for CSV.
- `source/statistics.R`, `covariance.R`, `utility_oracle.R` are byte-preserved from the accepted reproduction's `source/`. Read them with 06: observed statistic, covariance bound/decay, Oracle draws. The author's accepted evidence remains unchanged under `_local/oracle-table-figure/results/`; no historical-number adjudication is required to use this code.

The saved-input route is the accepted reproduction boundary. Packaging fixtures and saved-task summary checks are separate execution evidence, not a new scientific validation. The generic HiGHS runner and operational OSCAR scripts are retained only in ignored `_local/packaging-backup-20260930/`.
