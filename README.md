# Smartphone addiction, source-specific social support and anxiety

R code for the manuscript's descriptive statistics, measurement analyses, sex-grouped path models, GAMLSS and Bayesian longitudinal networks.

## Data access

The de-identified participant-level datasets are deposited in ScienceDB at [10.57760/sciencedb.psych.000yz](https://doi.org/10.57760/sciencedb.psych.000yz). Access is restricted; requests should be submitted through ScienceDB.

## Run analyses

1. Install the required R packages and CmdStan (see [software requirements](docs/software.md)).
2. Copy `config.R.example` to `local_config.R`. Place the controlled-access `cross_sectional.csv` and `longitudinal.csv` in `data_raw/`, or set their paths under `shared`. Set output locations as needed. The default input mode is `shared`.
3. From this repository, run:

```r
source("run_pipeline.R", encoding = "UTF-8")
run_pipeline("plan")
run_pipeline("all")
```

`all` prepares the datasets and runs the listed analysis stages to generate model results, tables, figures and diagnostics. See [manuscript outputs](docs/outputs.md) for the output map and software-specific limitations. Bayesian fitting and path-model bootstrap resampling can take substantial time. `plan` lists the stages and input files without fitting.

Stages can be run separately:

```r
run_pipeline("prepare")
run_pipeline("descriptive")
run_pipeline("measurement")
run_pipeline("binary")
run_pipeline("path_tables")
run_pipeline("gamlss")
run_pipeline("gaussian_fit")
run_pipeline("gaussian_results")
run_pipeline("gaussian_ppc")
run_pipeline("ordinal_fit")
run_pipeline("ordinal_results")
```

Participant-level data and fitted models are not included in this code repository. Input files remain in `data_raw/`; prepared datasets are written to `data/` and results to `results/`. Participant files, generated results and `local_config.R` are excluded from version control. See [input data and scoring](docs/data.md) for coding and linkage rules. For source questionnaire exports, set `input_mode = "raw"` and configure the `raw` file paths, worksheet names and column mappings. Excel and CSV are supported.

MCMC and bootstrap output can vary across software and hardware environments. Analysis runs record software versions, seeds and configuration.

## Resume and parallel execution

Gaussian and ordinal caches are checked against the current data, formulas, software versions and sampling specification. Compatible fits are reused. Incompatible caches stop execution; use a different results directory to keep both analyses.

Chain concurrency is configured with `network$cores`. Ordinal targets can also be divided between independent terminal processes:

```sh
Rscript run_pipeline.R ordinal_fit local_config.R support
Rscript run_pipeline.R ordinal_fit local_config.R mpiq
```

Target lists such as `kin,offline` are accepted. Use disjoint target lists for concurrent processes. Compatible cached targets are skipped. Run `ordinal_results` once after all targets finish. A target lock prevents two processes from writing the same fit. If a process is terminated, check that no fit is running before manually removing its empty `.lock` directory.

## Checks and outputs

```sh
Rscript run_code_checks.R
Rscript tests/test_metrics.R
Rscript tests/test_data_pipeline.R
```

Tests use synthetic questionnaire records and predefined network examples to check data preparation and network calculations without participant data. Data-preparation tests cover two sample sizes and write test files to `results/tests/`. The GitHub Actions workflow is configured to run source, network-metric and data-preparation checks on Linux. These checks do not fit Bayesian models. See [manuscript outputs](docs/outputs.md) for the manuscript figure/table map and generated files, and [analysis definitions](docs/methods.md) for model and metric definitions.
