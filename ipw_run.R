#!/usr/bin/env Rscript
####################################################################################################
# Description: This is the primary script for running the IPW analysis for Sepsis Time Zero
#
# DEPENDENT SCRIPTS
#   ipw_setup.R
#   ipw_boot.R
#
# INPUT FILES
#   $INPUT_FILE.{rds|csv}  <- the input person-time dataset
#
# OUTPUT FILES
#   $INPUT_FILE.after_setup.rds  <- a cache of the dataset after cleanup/setup procedures
#   $INPUT_FILE.adj-$MODEL_NAME.est_debug.txt  <- contains debug info for the point estimates
#   $INPUT_FILE.adj-$MODEL_NAME.est_summary.csv  <- a summary of point estimates
#   $INPUT_FILE.adj-$MODEL_NAME.boot-$NUM_BOOTSTRAPS.rds  <- estimates for the bootstrap replicates
#   $INPUT_FILE.adj-$MODEL_NAME.boot_summary-$NUM_BOOTSTRAPS.csv  <- a summary of the bootstraps
#
# Authors: Theodore Pak, Anna Agan, Jessica G. Young
####################################################################################################

require(tools)
require(data.table)
require(progressr)

source("ipw_setup.R")
source("ipw_boot.R")

usage <- "
USAGE
ipw_run.R MODEL_NAME INPUT_FILE [NUM_WORKERS] [NUM_BOOTSTRAPS]

ARGUMENTS
- MODEL_NAME is the adjustment model to use. Valid options are:
    - max
    - max_HCA
    - baseline
    - baseline_HCA
    - none
- INPUT_FILE is the path to the person-time dataset.
- NUM_WORKERS is the number of workers to run in parallel for computing
  bootstraps. If none is given, parallelization is disabled.
- NUM_BOOTSTRAPS is the number of bootstraps to run. If not given a default
  of 1000 is used.

OUTPUT
Output is saved in the same directory as the input file, with the same 
base name as the input file, but with different suffixes for all the various
outputs."

args <- commandArgs(trailingOnly = TRUE)

# Check arguments and stop with a usage message if they are invalid
if (length(args) < 2) {
    stop(paste0("Not enough arguments were given.\n", usage))
} else if (length(args) >= 2) {
    adj_model <- args[1]
    input_file <- args[2]
}
n_workers <- ifelse(length(args) >= 3, as.numeric(args[3]), 1)
n_reps <- ifelse(length(args) >= 4, as.numeric(args[4]), 1000)
if (length(args) > 4) {
    stop(paste0("A maximum of three arguments are allowed.\n", usage))
}

out_prefix <- file_path_sans_ext(input_file)
out_after_setup_rds <- paste0(out_prefix, ".after_setup.rds")

if (file.exists(out_after_setup_rds)) {
    message(paste0("Reading cached post-setup dataset from ", out_after_setup_rds))
    dt.pt_ipw <- readRDS(out_after_setup_rds)
} else {
    # Open and read the INPUT_FILE
    message("Reading input file ", input_file)
    if (toupper(file_ext(input_file)) == "RDS") {
        dt.pt <- readRDS(input_file)
    } else if (toupper(file_ext(input_file)) == "CSV") {
        dt.pt <- fread(input_file)
    }

    message("Setting up dataset.")
    dt.pt_ipw <- ipw_setup(dt.pt)
    message(paste0("Saving post-setup dataset to ", out_after_setup_rds))
    saveRDS(dt.pt_ipw, out_after_setup_rds)
    message("Saved.")
    rm(dt.pt)
}

# Constants used in the remainder of the code
Kint <- 12                # Number of time intervals
maxtime <- Kint - 1       # Max time interval is (Kint - 1) since time starts at 0
cutTimes <- c(0:maxtime)  # A vector with all valid time intervals
# A two-column matrix with start and end intervals for each time window to compare
compareTimes <- t(matrix(c(0, 5, 0:5, 2, 5, 6, 11, 6:11), nrow=2))

# Setup parallel processing, if enabled by the NUM_WORKERS option
set.seed(2000)
if (n_workers > 1) {
    require(doFuture)
    require(doRNG)
    registerDoFuture()
    registerDoRNG()
    if (supportsMulticore()) { 
        plan(multicore, workers = n_workers) 
    } else { 
        # R on windows doesn't allow for multicore (forking) and shared memory;
        # instead, use multisession with up to 16GB in globals exported to each worker
        plan(multisession, workers = n_workers)
        options(future.globals.maxSize = 16 * 1024^3)
    }
}
options(progressr.enable = TRUE)
handlers(handler_progress(format = "[:bar] :percent (:current/:total) eta: :eta"))

message("Creating point estimates.")
out_est_debug <- paste0(out_prefix, ".adj-", adj_model, ".est_debug.txt")
dt.est <- ipw_est(dt.pt_ipw, Kint, cutTimes, compareTimes, adj = adj_model,
    out_debug = out_est_debug)
out_est_summary <- paste0(out_prefix, ".adj-", adj_model, ".est_summary.csv")
dt.est.summary <- summarize_ipw_est(dt.est, compareTimes)
fwrite(dt.est.summary, out_est_summary)
message(paste0("Saved point estimate summary to ", out_est_summary))
message(paste0("Saved point estimate debug information to ", out_est_debug))

message(paste0("Running ", n_reps, " bootstraps."))
dt.boot.adj <- with_progress({
    ipw_boot(dt.pt_ipw, Kint, cutTimes, compareTimes, adj = adj_model, n_reps = n_reps,
        .parallel = n_workers > 1)
})

out_boot_raw <- paste0(out_prefix, ".adj-", adj_model, ".boot-", n_reps, ".rds")
saveRDS(dt.boot.adj, out_boot_raw)
message(paste0("Saved raw bootstrap output to ", out_boot_raw))

out_boot_summary <- paste0(out_prefix, ".adj-", adj_model, ".boot_summary-", 
    n_reps, ".csv")
fwrite(summarize_ipw_boot(dt.boot.adj, dt.est.summary, compareTimes), out_boot_summary)
message(paste0("Saved bootstrap summary statistics to ", out_boot_summary))
