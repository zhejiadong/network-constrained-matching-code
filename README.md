# Code for simulation studies and data application of FHS

The repository has two parts:

- `simulation/`: data generation, matching, null inference, sensitivity, power and matching-evaluation figures.
- `fhs_application/`: preparation of authorized FHS/dbGaP files followed by the application analysis.

Matching and inference use [netmatchRI](https://github.com/zhejiadong/netmatchRI).

## Install

```r
install.packages(c("remotes", "ggplot2", "mvtnorm", "ape"))
remotes::install_github("zhejiadong/netmatchRI", upgrade = "never")
```

The simulation uses Gurobi, as in the paper; its R package and a valid license
are required. To use HiGHS instead, change `solver = "gurobi"` to
`solver = "highs"` in `simulation/helpers.R`. The FHS analysis uses HiGHS.

## Simulation

Run these commands in `simulation/`:

```sh
Rscript 01_generate_data.R
Rscript 02_run_analysis.R
Rscript 03_make_tables_figures.R
```

This runs the null experiment with 300 units, three dependence settings and
500 replications. For the power analysis, add `power` to each command:

```sh
Rscript 01_generate_data.R power
Rscript 02_run_analysis.R power
Rscript 03_make_tables_figures.R power
```

Data are saved in `data/`, tables in `results/`, and figures in `figures/`.
Rerunning a step replaces its outputs.

| Output | Content |
| --- | --- |
| `null_table.csv`, `null_table.tex` | Null rejection rates |
| `cov_mahal.pdf`, `net_avg.pdf` | Matching diagnostics under strong dependence |
| `ri_grid_linePlot.pdf` | Sensitivity analysis |
| `matching_eval_all_combined.pdf` | Supplementary matching diagnostics |
| `power_figure.pdf` | Power analysis |

## FHS application

See [fhs_application/README.md](fhs_application/README.md) for data requirements
and run instructions. FHS data are not included.
