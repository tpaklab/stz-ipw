# Optimizing Sepsis Time Zero Definitions Based on Associations Between Time-to-Antibiotics and Mortality

This repository contains R code for the inverse probability weighting (IPW) analysis within the above manuscript.

## Prerequisites

R version ≥4.0.2 and the following packages: `data.table`, `progressr`, and `tidyverse`.

## File descriptions

### Scripts

1. `ipw_run.R` 

    **Description:** The primary script that runs the IPW analysis, including loading of the input person-time dataset, data cleanup, generation of point estimates, bootstrapping, generating summary tables, and saving output files. This script can be invoked directly and expects 2 to 4 command line arguments:

    ```
    $ `ipw_run.R` MODEL_NAME INPUT_FILE [NUM_WORKERS] [NUM_BOOTSTRAPS]
    ```

    `MODEL_NAME` specifies the variables included in the pooled over time logistic regression model. Possible values include `none`, `baseline_HCA`, and `max_HCA` (the fully adjusted model including time-varying covariates).

    **Input files:**

    - `$INPUT_FILE.{rds|csv}`  <- the input person-time dataset
    
    **Output files:**

    - `$INPUT_FILE.after_setup.rds`  <- a cache of the dataset after cleanup/setup procedures
    - `$INPUT_FILE.adj-$MODEL_NAME.est_debug.txt`  <- contains debug info for the point estimates
    - `$INPUT_FILE.adj-$MODEL_NAME.est_summary.csv`  <- a summary of point estimates
    - `$INPUT_FILE.adj-$MODEL_NAME.boot-$NUM_BOOTSTRAPS.rds`  <- estimates for the bootstrap replicates
    - `$INPUT_FILE.adj-$MODEL_NAME.boot_summary-$NUM_BOOTSTRAPS.csv`  <- a summary of the bootstraps

### Library files

1. `ipw_setup.R`

    **Description:** Setup and cleanup functions for the input person-time datasets for ipw_run.R. The most important exported function is `ipw_setup()`, which runs the entire setup process.

2. `ipw_boot.R`

    **Description:** Contains the primary steps of the IPW analysis that are called by `ipw_run.R`.

    The most important exported functions are: 

    - `ipw_est()`  <- creates point estimates 
    - `ipw_boot()`   <- creates bootstrapped estimates that can be used to construct a 95% CI
    - `summarize_ipw_est()`  <- creates summary statistics for the point estimates
    - `summarize_ipw_boot()`  <- creates summary statistics for the bootrapped estimates

## License

MIT. See `LICENSE.txt`

## Contact 

For questions, clarifications, or further support, please contact:

Theodore Pak, MD, PhD — theodore (dot) pak (at) uci (dot) edu