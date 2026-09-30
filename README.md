# Network-dependent causal inference: paper code

Original paper code for the simulation study and Framingham Heart Study (FHS) application, organized for line-by-line author review.

The repository has two independent parts:

- `simulation/`: original network/DGP, Gurobi matching, null inference (including Oracle), sensitivity, parameter selection, power and matching-evaluation figures.
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

Alternatively, accepted saved tasks may be placed under ignored `simulation/tasks/` and summarized without recomputation. See `simulation/README.md` for the complete 01–13 read order, including the original power pipeline and matching evaluation. Locally generated simulation CSVs are permitted but must remain ignored.

## FHS application

The FHS data are controlled-access and are not distributed. Authorized users place the required dbGaP files under `fhs_application/data_raw/`, then run from `fhs_application/`:

```text
Rscript 01_prepare_data.R
Rscript 02_analysis.R
```

`01_prepare_data.R` specifies tables/variables and creates a local analysis object. `02_analysis.R` uses it for matching, inference, sensitivity and plots. Keep all controlled inputs and outputs inside the authorized environment.

## Dependencies

Checked using R 4.5.0 and existing libraries; nothing was installed. Direct dependencies and observed versions:

- Saved-input inference: base/recommended R only; plotting: `ggplot2` 4.0.3. The standalone `draw_Z_oracle()` helper additionally needs `mvtnorm` 1.3.3 when called without precomputed Gaussian draws; the runner supplies those draws.
- Original network/DGP scripts: `igraph` 2.1.4, `tidyr` 1.3.1, `data.table` 1.17.8, `dplyr` 1.1.4, `energy` 1.7.12, `sna` 2.8, `mvtnorm` 1.3.3, `RItools` 0.3.5, `ape` 5.8.1, `Matrix` 1.7.3. Historical library calls are retained, including unused imports.
- Original matching: `dplyr`, `Matrix`, `igraph`, `MASS` 7.3.65, `gurobi` 10.0.1 (separate installation and valid solver license), `R.utils` 2.13.0, `optmatch` 0.10.8, `foreach` 1.5.2, `doParallel` 1.0.17. Matching requires the solver; saved-input inference does not.
- FHS: `igraph`, `ape`, `netmatchRI` 0.1.0, `ggplot2`, `MASS`, `highs` 1.14.0.2, plus those packages' declared dependencies. FHS dual matching explicitly selects HiGHS.
- Power/selection: `dplyr`, `foreach`, `doParallel`; evaluation/plots additionally use `patchwork`, `scales`, `tidyr`, `MASS`. `patchwork` availability was checked.

[`netmatchRI`](https://github.com/zhejiadong/netmatchRI) remains the public reusable method package. The FHS analysis calls its API. The simulation deliberately retains the exact legacy matching and accepted inference kernels rather than substituting package calls: this is a paper-reproduction snapshot, not a fork or proposed package change. Versions above describe the packaging environment, not a historical lockfile.

## Verification boundary

Accepted original-input OSCAR null validation is reused, not rerun. Packaging checks parse every R file, preserve source hashes, export saved null grid rows and execute the original selector. Power plots and one original power helper were exercised only on artificial fixtures. Full power, generation/matching, study-data matching evaluation and controlled FHS execution were not rerun; historical bit-identical results are not promised. All publication-bound files require author review and a license decision.
