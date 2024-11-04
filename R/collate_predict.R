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
source("R/functions_predict.R")
config <- config::get(config = config_name)

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

read_path <- get_path("read", config_name, task_id)

density_tasks <- list.files(read_path)
message("Tasks to collate ", length(density_tasks))


pb <- txtProgressBar(max = length(density_tasks), style = 1)
# for(i in seq_along(density_tasks)){
for(i in 1:1){

  task_id <- density_tasks[i]

  rds_file <- file.path(read_path, task_id, "simulation_data.rds")

  if(file.exists(rds_file)){
    rds <- read_rds(rds_file)
  } else {
    next
  }


  psrf <- rds$psrf |>
    as_tibble() |>
    mutate(node_names = rownames(rds$psrf)) |>
    filter(node_names != "psi_phi")

  bad_mcmc <- rds$bad_mcmc | any(psrf$`Upper C.I.` > 1.1)

  if(bad_mcmc) next

  start_density <- rds$start_density

  samples <- rds$posterior_samples |>
    select(
      contains("beta1"),
      contains("beta_p"),
      contains("log_gamma["),
      contains("log_rho["),
      contains("p_mu["),
      contains("log_nu"),
      contains("phi_mu"),
      contains("psi_phi")
    )

  constants <- rds$constants
  data <- rds$data

  ls <- data_posteriors(samples, constants, data)

  print(str(ls))

  # |>
  #   add_ids(t_id, dens)

  setTxtProgressBar(pb, i)
}
close(pb)










