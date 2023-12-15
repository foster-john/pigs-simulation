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

library(parallel)
library(nimble)
library(coda)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)

config_name <- "dev"
config <- config::get(config = config_name)

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

args <- commandArgs(trailingOnly = TRUE)
n_threads <- args[1]
message("Number of cores available ", n_threads)

source("R/run_simulation.R")
cl <- makeCluster(as.numeric(n_threads))
sim <- run_simulation(cl, config, df)
stopCluster(cl)

# -----------------------------------------------------------------
# Summarize output ----
# -----------------------------------------------------------------

collate_mcmc_output(config, sim)



