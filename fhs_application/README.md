# Framingham Heart Study application

The analysis uses Exam 4 obesity, Exam 5 height-indexed left ventricular mass,
and baseline covariates from the FHS Offspring cohort.

## Data

Place the authorized dbGaP files in `data_raw/`. Use both consent groups
(`c1` and `c2`) with the original `.txt.gz` filenames.

| Table | Filename prefix before `c1.` or `c2.` |
| --- | --- |
| Exam 4 clinic | `phs000007.v35.pht000033.v13.p16.` |
| Exam 4 echocardiography | `phs000007.v35.pht000300.v12.p16.` |
| Exam 5 clinic | `phs000007.v35.pht000034.v12.p16.` |
| Exam 5 echocardiography | `phs000007.v35.pht000301.v12.p16.` |
| Dates and demographics | `phs000007.v35.pht003099.v10.p16.` |
| Social network | `phs000153.v9.pht000836.v8.p8.` |

## Run

Install the packages listed in the repository README, then run in `fhs_application/`:

```sh
Rscript 01_prepare_data.R
Rscript 02_analysis.R
```

The first script prepares the cohort and network. The second runs matching,
inference and sensitivity analysis with `kappa <- 3`. For the supplementary
analyses, change `kappa` to 2 or 4 at the top of `02_analysis.R`, then rerun it.

Prepared data are saved in `data/`, numerical results in `results/`, and
sensitivity, balance and covariate-similarity plots in `figures/`.
Output filenames identify the kappa value.

FHS data require authorized access and are not distributed here. Keep all
participant-level inputs and outputs in the authorized environment.
