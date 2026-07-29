# Smartphone addiction and anxiety in Chinese university students

This repository contains the analysis code associated with a study of
smartphone addiction, anxiety, source-specific social support and
sex-specific differences among Chinese university students.

## Repository status

This is currently a private working repository. The analysis code and
documentation are being organized and checked in preparation for a
reproducible research release.

## Analyses

The repository will contain code for:

- data preparation and variable construction
- descriptive statistics and group comparisons
- mediation and moderation analyses
- generalized additive models for location, scale and shape (GAMLSS)
- Bayesian cross-lagged panel network analyses
- posterior network comparisons
- main and supplementary tables and figures
- sensitivity and robustness analyses

## Data availability

Participant-level data from the CPMC dataset are not included in this
repository because access is restricted by the consortium's data-use,
research ethics and participant privacy requirements.

Researchers who wish to use the data must follow the formal CPMC data
application and approval procedures.

## Preregistration

The hypotheses and planned analyses were preregistered on Open Science
for Psychology.

- preRegisterID: 202511.00033

The preregistration was completed after data collection because the CPMC
dataset had been collected independently of the present study. It was
completed before the reported analyses were conducted.

## Authors and maintainers

- Li-Kun Ge — first author and repository maintainer
- Gao-Xia Wei — corresponding author
- Additional contributors will be listed according to their contributions

## Planned repository structure

```text
R/
├─ 00_setup.R
├─ 01_data_preparation.R
├─ 02_descriptive_statistics.R
├─ 03_mediation_moderation.R
├─ 04_gamlss_analysis.R
├─ 05_gamlss_figures.R
├─ 06_bayesian_network_analysis.R
├─ 07_network_comparisons.R
├─ 08_main_outputs.R
└─ 09_supplementary_outputs.R

docs/
synthetic_data/
outputs/
```

## Reproducibility

Participant-level data are not publicly distributed. The repository will
provide:

- complete analysis scripts
- software and package requirements
- variable descriptions
- an analysis workflow
- a manuscript output-to-code map
- synthetic or simulated data where appropriate for testing the workflow