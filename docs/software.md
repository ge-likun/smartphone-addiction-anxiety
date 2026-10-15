# Software

The reference software environment comprises R 4.5.1, brms 2.22.0, cmdstanr 0.9.0, CmdStan 2.36.0, gamlss 5.4-22 and gamlss.dist 6.1-1. Use these versions for a comparable software environment. Scripts record the versions actually installed; caches from different environments are not silently reused.

Required packages:

```r
install.packages(c("readxl", "dplyr", "tibble", "tidyr", "lavaan", "gamlss",
                  "gamlss.dist", "brms", "posterior", "rstan", "bayesplot",
                  "qgraph", "ggplot2", "patchwork", "cowplot", "officer",
                  "flextable", "writexl", "svglite", "ragg", "jsonlite", "magick"))
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev",
                                       getOption("repos")))
cmdstanr::install_cmdstan(version = "2.36.0")
```

On Windows, CmdStan requires a compatible C++ toolchain. Set `CMDSTAN` to an existing installation if it is not found by cmdstanr.

Automatic package installation is not performed by analysis scripts. Output directories record `sessionInfo()`, seeds and configuration information.
