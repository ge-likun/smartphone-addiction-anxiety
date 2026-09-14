# Smartphone addiction, social support and anxiety

This repository contains the R code used for three parts of the manuscript:
descriptive and sex-grouped path analyses, GAMLSS distributional regression,
and Bayesian longitudinal network analysis.

## Software

The reference analyses used R 4.5.1. The descriptive and path-model workflow
used `dplyr` 1.1.4, `lavaan` 0.6-21, `flextable` 0.9.7 and `officer` 0.6.7.
The GAMLSS workflow used `gamlss` 5.4-22 and `gamlss.dist` 6.1-1. The Bayesian
network model used `brms` 2.22.0, `cmdstanr` 0.9.0 and CmdStan 2.36.0. Required
packages are checked by the analysis scripts and are not installed
automatically.

## Repository structure

### `01_descriptive_and_path_models`

| File | Purpose |
|---|---|
| `00_functions.R` | Shared functions for effect sizes, equivalence tests and publication tables |
| `01_descriptive_tables.R` | Main Table 1 and Supplementary Tables 1 and 2 |
| `02_moderated_mediation.R` | Sex-grouped conditional-association models and Supplementary Table 5 |
| `03_retest_stability.R` | Retest-stability summaries |
| `run_all.R` | Runs the four scripts in order |

### `02_gamlss`

| File | Purpose |
|---|---|
| `00_zanbi_order_safe.R` | Order-safe ZANBI functions and preflight checks for `gamlss.dist` 6.1-1 |
| `01_setup_and_data.R` | Package checks, data validation, analysis strata and model formulas |
| `02_fit_models.R` | Fits and verifies NO, ZAGA, ZAIG, ZANBI and ZABB models |
| `03_residual_diagnostics.R` | Reproducible normalized quantile residuals and boundary summaries |
| `04_figures.R` | Main Figures 2 and 3 and Supplementary Figures 4–6 |
| `05_tables.R` | Supplementary Tables 6–8 |
| `06_covariate_sensitivity.R` | Supplementary Figure 7 |
| `run_all.R` | Runs the seven scripts in order |

The GAMLSS code combines the verified order-safe ZANBI implementation used for
the reported discrete sensitivity results with data-fingerprinted model cache
names. Information criteria are compared only within the continuous/mixed set
(NO, ZAGA and ZAIG) and within the discrete set (ZANBI and ZABB).

### `03_bayesian_network`

| File | Purpose |
|---|---|
| `00_packages.R` | Checks the required R packages |
| `01_config.R` | Defines nodes, groups, seeds, paths and sampling settings |
| `02_data.R` | Imports, recodes and validates the two-wave network data |
| `03_model.R` | Defines, fits and loads the multivariate `brms` model |
| `04_diagnostics.R` | Summarizes convergence and sampler diagnostics |
| `05_networks.R` | Reconstructs group-specific posterior network matrices |
| `06_metrics.R` | Calculates network summaries, contrasts and distances |
| `07_figures.R` | Creates Main Figures 4 and 5 |
| `08_tables.R` | Creates Supplementary Tables 9–12 |
| `run_fit_model.R` | Fits a new Bayesian model |
| `run_postprocessing.R` | Generates results from an existing fitted model |

The separate Bayesian scripts correspond to distinct stages of one analysis.
This keeps the computationally expensive model fit independent of diagnostics,
network reconstruction, figures and tables.

## Input files

Participant-level data and fitted models are not included. The default input
files are `analysis_data.csv`, `retest_data.csv` and `network_data.csv` in the
repository root. These filenames are excluded by `.gitignore`.

Alternative locations can be supplied before running the code:

```r
Sys.setenv(
  ANALYSIS_DATA_FILE = "path/to/analysis_data.csv",
  RETEST_DATA_FILE = "path/to/retest_data.csv",
  GAMLSS_DATA_FILE = "path/to/analysis_data.csv",
  NETWORK_DATA_FILE = "path/to/network_data.csv",
  NETWORK_MODEL_FILE = "path/to/brms_model.rds"
)
```

The scripts calculate sample sizes from the supplied data. Study-design checks
for score ranges, grouping variables, analysis strata and sufficient variation
are retained.

## Running the code

Set the working directory to the relevant analysis folder and source its runner.

```r
# Descriptive, path-model and retest analyses
setwd("01_descriptive_and_path_models")
source("run_all.R", encoding = "UTF-8")

# GAMLSS analyses
setwd("../02_gamlss")
source("run_all.R", encoding = "UTF-8")

# Bayesian post-processing using an existing fitted model
setwd("../03_bayesian_network")
source("run_postprocessing.R", encoding = "UTF-8")
```

The Bayesian model is fitted only when `run_fit_model.R` is sourced. Model
fitting is computationally intensive. The post-processing scripts use an
existing fitted model and do not start a new fit.

Generated files are written under `results/`, which is excluded from version
control. Each workflow writes software information with its results.
