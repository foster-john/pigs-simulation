# --------------------------------------------------------------------
#
# Script for collating predicted values
#
# John Foster
#
# --------------------------------------------------------------------

message("\n=== PARAMETERS ===\n")

# will get config_name from here
source("R/functions_collate.R")
config <- config::get(config = config_name)

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

read_path <- get_path("read", config_name, task_id)

density_tasks <- list.files(read_path)
message("Tasks to collate ", length(density_tasks))

nodes <- "parameters"
tasks_ls <- get_tasks(density_tasks, read_path, nodes)
all_samples <- tasks_ls$all_sample

all_param_samples <- all_samples |>
  select(
    simulation,
    start_density,
    contains("beta1"),
    contains("beta_p"),
    contains("log_gamma["),
    contains("log_rho["),
    contains("p_mu["),
    contains("log_nu"),
    contains("phi_mu"),
    contains("psi_phi")
  )

print(glimpse(all_param_samples))


