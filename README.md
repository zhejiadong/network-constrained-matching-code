# Network-dependent causal inference: paper code

Code-only author-review candidate for the simulation study and Framingham Heart Study (FHS) application. The author has accepted the completed OSCAR original-input null reproduction. This packaging preserves that lineage; it does not replace manuscript values or claim identical historical Monte Carlo streams.

The repository has two independent parts:

- `simulation/`: original network/DGP and Gurobi matching entrypoints, accepted original-input null inference (including Oracle), and rejection-rate summaries/figure.
- `fhs_application/`: preparation of authorized FHS/dbGaP files followed by the application analysis.

No data, CSV files, cached R objects, binaries, results or plots are distributed. The deny-by-default `.gitignore` exposes only the reviewed code and documentation. `_local/` is private history/evidence, not a publication input. License: pending author choice; no license grant is made here.

## Simulation

Read `simulation/README.md` for exact execution order and input schema. From `simulation/`, accepted saved synthetic `inputs/null_dep1.rds` through `null_dep3.rds` support:

```text
Rscript 06_run_null_inference.R 1 smoke
```

The `smoke` argument processes one dataset with the unchanged 500 Oracle draws. For a full null inference run, execute tasks 1 through 150 without `smoke`, then:

```text
Rscript 07_summarize.R
```

Alternatively, previously accepted task outputs may be placed under ignored `simulation/tasks/` and summarized without recomputation. Scripts 01–05 include the original generator/matching lineage for new null datasets; this packaged end-to-end route has not been run in full. New solver outputs are not asserted to equal frozen historical matches. No power-reproduction claim is made.

## FHS application

The FHS data are controlled-access and are not distributed. Authorized users place the required dbGaP files under `fhs_application/data_raw/`, then run from `fhs_application/`:

```text
Rscript 01_prepare_data.R
Rscript 02_analysis.R
```

`01_prepare_data.R` specifies tables/variables and creates a local analysis object. `02_analysis.R` uses it for matching, inference, sensitivity and plots. These scripts were parsed, not executed against controlled inputs in this packaging task. A synthetic tab-delimited intake fixture was tested; that is not verification of FHS cohort construction or numerical results. Keep all FHS inputs, IDs, edges, matrices, matches, results, figures and console logs inside the authorized environment, never in a public repository.

## Dependencies

Checked using R 4.5.0 and existing libraries; nothing was installed. Direct dependencies and observed versions:

- Saved-input inference: base/recommended R only; plotting: `ggplot2` 4.0.3. The standalone `draw_Z_oracle()` helper additionally needs `mvtnorm` 1.3.3 when called without precomputed Gaussian draws; the runner supplies those draws.
- Original network/DGP scripts: `igraph` 2.1.4, `tidyr` 1.3.1, `data.table` 1.17.8, `dplyr` 1.1.4, `energy` 1.7.12, `sna` 2.8, `mvtnorm` 1.3.3, `RItools` 0.3.5, `ape` 5.8.1, `Matrix` 1.7.3. Historical library calls are retained, including unused imports.
- Original matching: `dplyr`, `Matrix`, `igraph`, `MASS` 7.3.65, `gurobi` 10.0.1 (separate installation and valid solver license), `R.utils` 2.13.0, `optmatch` 0.10.8, `foreach` 1.5.2, `doParallel` 1.0.17. Matching requires the solver; saved-input inference does not.
- FHS: `igraph`, `ape`, `netmatchRI` 0.1.0, `ggplot2`, `MASS`, `highs` 1.14.0.2, plus those packages' declared dependencies. FHS dual matching explicitly selects HiGHS.

[`netmatchRI`](https://github.com/zhejiadong/netmatchRI) remains the public reusable method package. The FHS analysis calls its API. The simulation deliberately retains the exact legacy matching and accepted inference kernels rather than substituting package calls: this is a paper-reproduction snapshot, not a fork or proposed package change. Versions above describe the packaging environment, not a historical lockfile.

## Verification boundary

Completed OSCAR validation is reused, not rerun: original matched null inputs, three dependence settings, 500 datasets each, 500 seeded Oracle draws per design, original 5-by-8 sensitivity grid. Packaging tests separately parse every R file, exercise deterministic artificial fixtures and checkpoint reuse, and summarize saved accepted tasks without new simulation. Full FHS execution and a new end-to-end generator/matching run are not claimed verified. All publication-bound files still require the author's line-by-line review.
