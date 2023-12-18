# --------------------------------------------------------------------
#
# Workflow script for pig removal simulations
#
# John Foster
#
# General workflow (conducted with run_simulation):
#
# 1. Load config and MIS data
#
# 2. run_simulation
#   - Simulate/bootstrap data
#     - 1-method properties
#     - 2- to 5-method properties
#   - Simulate eco/take dynamics
#   - Fit MCMC
#   - Check MCMC
# 3. Summarize output
#
# --------------------------------------------------------------------

# renv::load("/home/john.foster/pigs-simulation/")
renv::config(renv.config.sandbox.enabled = FALSE)

config_name <- "hpc_production"
config <- config::get(config = config_name)

library(nimble)
library(parallel)
library(coda)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)


# -----------------------------------------------------------------
# Load MIS data ----
# -----------------------------------------------------------------
message(" MIS data intake")
top_dir <- config$top_dir
data_dir <- config$data_dir
df <- read_rds(file.path(top_dir, data_dir, "insitu/MIS_4weekPP.rds"))
df <- df |>
  select(-method) |>
  rename(property = agrp_prp_id,
         method = Method)

# -----------------------------------------------------------------
# Run simulation ----
# -----------------------------------------------------------------

message("Number of cores available ", n_threads)
source("R/run_simulation.R")

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
run_simulation(config, df, task_id)

message("\n\nDONE!")


