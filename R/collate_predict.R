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

  samples <- rds$posterior_samples
  constants <- rds$constants
  constants$samples <- as.matrix(samples)
  data <- rds$data

  message("NIMBLE")
  print(str(constants))
  print(str(data))

  ls <- data_posteriors(samples, constants, data)

  N <- rds$N
  take <- rds$take


  message("Take")

  print(str(ls))
  print(str(take))

  y_pred <- t(as.matrix(ls$y_pred))
  colnames(y_pred) <- paste0("iter_", 1:ncol(y_pred))

  print(str(y_pred))
  
  y_pred <- y_pred |>
    as_tibble() |>
    mutate(PPNum = take$PPNum,
           N = take$N,
           take = take$take,
           property = take$property)

  left_join(y_pred, take)



  # |>
  #   add_ids(t_id, dens)

  setTxtProgressBar(pb, i)
}
close(pb)










