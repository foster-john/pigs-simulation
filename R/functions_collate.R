library(stringr)
library(dplyr)
library(tidyr)
library(nimble)
library(readr)
library(coda)

config_name <- "hpc_production"

get_path <- function(type, config_name, task_id){

  config <- config::get(config = config_name)
  top_dir <- config$top_dir
  out_dir <- config$out_dir
  analysis_dir <- config$analysis_dir
  dev_dir <- config$dev_dir
  model_dir <- config$model_dir
  project_dir <- config$project_dir

  start_vec <-
    c("0.3",
      "1.475",
      "2.65",
      "3.825",
      "5")
  start_density <- start_vec[task_id]
  density_dir <- paste0("density_", start_density)

  if(type == "read") path <- file.path(top_dir, project_dir, out_dir, dev_dir, model_dir, density_dir)
  if(type == "write") path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, model_dir, density_dir)
  return(path)
}

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])




add_ids <- function(df, task_id, s_density){
  df |>
    mutate(simulation = task_id,
           start_density = s_density)
}

bind_samples <- function(all_bind, ls, t_id, dens){
  samples <- ls$posterior_samples |>
    add_ids(t_id, dens)
  bind_rows(all_bind, samples)
}

bind_y <- function(all_bind, ls, t_id, dens){
  y_pred <- ls$posterior_take
  colnames(y_pred) <- 1:ncol(y_pred)
  y_pred <- y_pred |>
    as_tibble() |>
    add_ids(t_id, dens)
  bind_rows(all_bind, y_pred)
}

bind_take <- function(all_bind, ls, t_id, dens){
  take <- ls$take |>
    add_ids(t_id, dens) |>
    mutate(p_id = 1:n())
  bind_rows(all_bind, take)
}

bind_post_summaries <- function(all_bind, node, ls, t_id, dens){
  df <- ls[[node]]
  tb <- as_tibble(df) |>
    add_ids(t_id, dens) |>
    mutate(p_id = 1:n(),
           parameter = node)
  bind_rows(all_bind, tb)
}

bind_N <- function(all_bind, ls, t_id, dens){
  obs_flag <- ls$take |>
    select(property, county, PPNum) |>
    distinct() |>
    mutate(obs_flag = 1)

  N <- ls$N |>
    add_ids(t_id, dens) |>
    left_join(obs_flag) |>
    mutate(obs_flag = if_else(is.na(obs_flag), 0, obs_flag))
  bind_rows(all_bind, N)
}

bind_beta_p <- function(all_bind, ls, t_id, dens){
  # need a lookup table for known data model covariates
  bH <- tibble(
    method_idx = rep(1:nrow(ls$beta_p), ncol(ls$beta_p)),
    position = rep(1:ncol(ls$beta_p), each = nrow(ls$beta_p)),
    actual = as.numeric(ls$beta_p)
  ) |>
    add_ids(t_id, dens)
  bind_rows(all_bind, bH)
}

bind_methods <- function(all_bind, ls, t_id, dens){
  method_lookup <- ls$method_lookup |>
    add_ids(t_id, dens)
  bind_rows(all_bind, method_lookup)
}

bind_psrf <- function(all_bind, ls, t_id, dens){
  names <- rownames(ls$psrf)
  psrf <- ls$psrf |>
    as_tibble() |>
    mutate(node = names) |>
    add_ids(t_id, dens)
  bind_rows(all_bind, psrf)
}

select_pivot_longer <- function(df, node){
  df |>
    select(contains(node), simulation, start_density) |>
    pivot_longer(cols = -c(simulation, start_density),
                 names_to = "node")
}

recovered <- function(df){
  df |> mutate(parameter_recovered = if_else(actual >= low & actual <= high, 1, 0))
}

my_summary <- function(df){
  df |>
    summarise(low = quantile(value, 0.025),
              q1 = quantile(value, 0.25),
              med = quantile(value, 0.5),
              q3 = quantile(value, 0.75),
              high = quantile(value, 0.975),
              mu = mean(value),
              sd = sd(value))
}

recov_beta1 <- function(df, psrf){
  df |>
    group_by(simulation, node, method_idx, position, start_density) |>
    my_summary() |>
    left_join(all_beta_p) |>
    ungroup() |>
    recovered() |>
    left_join(psrf)
}

resid_beta1 <- function(df){
  df |>
    left_join(all_beta_p) |>
    mutate(value = value - actual) |>
    group_by(node, position, method_idx, start_density) |>
    my_summary() |>
    ungroup()
}

recov_beta_p <- function(df, psrf){
  df |>
    group_by(simulation, node, method_idx, position, start_density) |>
    my_summary() |>
    left_join(all_beta_p) |>
    ungroup() |>
    recovered() |>
    left_join(psrf)
}

resid_beta_p <- function(df){
  df |>
    left_join(all_beta_p) |>
    mutate(value = value - actual) |>
    group_by(node, position, method_idx, start_density) |>
    my_summary() |>
    ungroup()
}

recov_gamma <- function(df, H, psrf){
  df |>
    group_by(simulation, node, idx, start_density) |>
    my_summary() |>
    left_join(H) |>
    ungroup() |>
    recovered() |>
    left_join(psrf)
}

resid_gamma <- function(df, H){
  df |>
    left_join(H)|>
    mutate(value = value - actual) |>
    group_by(node, idx, start_density) |>
    my_summary() |>
    ungroup()
}

recov_phi <- function(df, known, psrf){
  df |>
    group_by(simulation, node, start_density) |>
    my_summary() |>
    mutate(actual = known) |>
    ungroup() |>
    recovered() |>
    left_join(psrf)
}

resid_phi <- function(df, known){
  phi_residual <- df |>
    mutate(actual = known) |>
    mutate(value = value - actual) |>
    group_by(node, start_density) |>
    my_summary() |>
    ungroup()
}

get_xn <- function(df, H){
  df |>
    select_pivot_longer("N[") |>
    filter(!is.na(value)) |>
    mutate(n_id = as.numeric(str_extract(node, "(?<=\\[)\\d*"))) |>
    left_join(H) |>
    filter(!is.na(abundance)) |>
    mutate(estimated_density = value / property_area)
}

get_post_take <- function(df, H){
  df |>
    pivot_longer(cols = -c(simulation, start_density),
                 names_to = "p_id") |>
    filter(!is.na(value)) |>
    mutate(p_id = as.numeric(p_id)) |>
    left_join(H)
}

take_calc <- function(df){
  df |>
    summarise(mae = mean(abs(value - take)),
              mpe = mean(abs((value+1) - (take+1))/(take+1))*100,
              mbias = mean(value - take),
              mse = mean((value - take)^2),
              rmse = sqrt(mse),
              nm_rmse = rmse / mean(take),
              sd_ratio = sd(value) / sd(take))
}

get_tasks <- function(density_tasks, path, nodes){

  if("all_samples" %in% nodes) all_samples <- tibble()
  if("all_take" %in% nodes)    all_take <- tibble()
  if("all_N" %in% nodes)       all_N <- tibble()
  if("all_beta_p" %in% nodes)  all_beta_p <- tibble()
  if("all_methods" %in% nodes) all_methods <- tibble()
  if("all_y" %in% nodes)       all_y <- tibble()
  if("all_area" %in% nodes)    all_area <- tibble()
  if("all_theta" %in% nodes)   all_theta <- tibble()
  if("all_p" %in% nodes)       all_p <- tibble()
  if("all_psrf" %in% nodes)    all_psrf <- tibble()

  message("Loop through tasks...")

  pb <- txtProgressBar(max = length(density_tasks), style = 1)
  for(i in seq_along(density_tasks)){

    task_id <- density_tasks[i]

    rds_file <- file.path(path, task_id, "simulation_data.rds")

    if(file.exists(rds_file)){
      rds <- read_rds(rds_file)
    } else {
      next
    }

    bad_mcmc <- rds$bad_mcmc | any(rds$psrf[,1] >= 1.1)
    converged <- rds$converged

    if(bad_mcmc) next

    start_density <- rds$start_density

    if("all_samples" %in% nodes) {
      all_samples <- bind_samples(all_samples, rds, task_id, start_density)
    }

    if("all_y" %in% nodes) {
      all_y <- bind_y(all_y, rds, task_id, start_density)
    }

    if("all_area" %in% nodes) {
      all_area <- bind_post_summaries(all_area, "posterior_potential_area", rds, task_id, start_density)
    }

    if("all_theta" %in% nodes) {
      all_theta <- bind_post_summaries(all_theta, "posterior_theta", rds, task_id, start_density)
    }

    if("all_p" %in% nodes) {
      all_p <- bind_post_summaries(all_p, "posterior_p", rds, task_id, start_density)
    }

    if("all_take" %in% nodes) {
      all_take <- bind_take(all_take, rds, task_id, start_density)
    }

    if("all_N" %in% nodes) {
      all_N <- bind_N(all_N, rds, task_id, start_density)
    }

    if("all_beta_p" %in% nodes) {
      all_beta_p <- bind_beta_p(all_beta_p, rds, task_id, start_density)
    }

    if("all_methods" %in% nodes) {
      all_methods <- bind_methods(all_methods, rds, task_id, start_density)
    }

    if("all_psrf" %in% nodes) {
      all_psrf <- bind_psrf(all_psrf, rds, task_id, start_density)
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  ls <- list()
  if("all_samples" %in% nodes) ls$all_samples <- all_samples
  if("all_take" %in% nodes)    ls$all_take <- all_take
  if("all_N" %in% nodes)       ls$all_N <- all_N
  if("all_beta_p" %in% nodes)  ls$all_beta_p <- all_beta_p
  if("all_methods" %in% nodes) ls$all_methods <- all_methods
  if("all_y" %in% nodes)       ls$all_y <- all_y
  if("all_area" %in% nodes)    ls$all_area <- all_area
  if("all_theta" %in% nodes)   ls$all_theta <- all_theta
  if("all_p" %in% nodes)       ls$all_p <- all_p
  if("all_psrf" %in% nodes)    ls$all_psrf <- all_psrf

  return(ls)

}


