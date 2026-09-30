# Framingham Heart Study application

This folder follows a code-only controlled-data workflow similar in disclosure level to Youjin Lee's public FHS repository: the required dbGaP/FHS tables, field names, eligibility rules, and analysis logic are documented, but no controlled data or participant-level derived objects are distributed.

## Use

1. Obtain authorized access to the required FHS/dbGaP release.
2. Place the compressed text tables listed in `01_prepare_data.R` under `data_raw/`.
3. Run:

```text
Rscript 01_prepare_data.R
Rscript 02_analysis.R
```

`01_prepare_data.R` is intended to create `data/fhs_analysis_data.rds` containing the eligible cohort, E4 network, and distance matrix. `02_analysis.R` runs covariate-only and dual-penalty matching, randomization inference, sensitivity analysis, and four application plots. Direct dependencies are listed in the root README.

These are the existing code-only application scripts, not a newly validated FHS reproduction. Packaging parsed both scripts and tested the tab-delimited intake helper with artificial rows only. Three literal-backslash-t delimiters in the header/unique-ID reader were corrected to actual tab delimiters; no cohort/analysis specification changed. Controlled-data execution, outputs and correspondence to manuscript figures remain unverified by this packaging task. Do not infer full FHS verification from a parse or synthetic intake test.

The `data_raw/`, `data/`, `results/`, and `figures/` directories are ignored. They may contain protected participant rows, identifiers, network edges, distance matrices, matched designs, or derived outcomes and must not be committed or transferred to a public repository.
