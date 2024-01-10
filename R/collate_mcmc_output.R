# --------------------------------------------------------------------
#
# Script for collating simulation results
#
# John Foster
#
# --------------------------------------------------------------------

library(stringr)
library(dplyr)
library(tidyr)
library(nimble)
library(readr)
library(coda)


config_name <- "hpc_production"
config <- config::get(config = config_name)

start_vec <-
  c("0.3",
    "1.475",
    "2.65",
    "3.825",
    "5")

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

start_density <- start_vec[task_id]

top_dir <- config$top_dir
out_dir <- config$out_dir
analysis_dir <- config$analysis_dir
dev_dir <- config$dev_dir
model_dir <- config$model_dir
project_dir <- config$project_dir

density_dir <- paste0("density_", start_density)

path <- file.path(top_dir, project_dir, out_dir, dev_dir, model_dir, density_dir)
# path <- "out/hpc/density_5"

density_tasks <- list.files(path)

message("Tasks to collate ", length(density_tasks))

# density_tasks <- 1:20
print(density_tasks)

all_samples <- tibble()
all_take <- tibble()
all_N <- tibble()
all_beta_p <- tibble()
all_methods <- tibble()
all_y <- tibble()
all_area <- tibble()
all_theta <- tibble()
all_p <- tibble()
all_psrf <- tibble()

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

  bad_mcmc <- rds$bad_mcmc | any(rds$psrf[,1] >= 1.2)
  converged <- rds$converged

  if(bad_mcmc) next

  start_density <- rds$start_density

  all_samples <- bind_samples(all_samples, rds, task_id, start_density)
  all_y <- bind_y(all_y, rds, task_id, start_density)
  all_area <- bind_post_summaries(all_area, "posterior_potential_area", rds, task_id, start_density)
  # all_theta <- bind_post_summaries(all_theta, "posterior_theta", rds, task_id, start_density)
  all_p <- bind_post_summaries(all_p, "posterior_p", rds, task_id, start_density)
  all_take <- bind_take(all_take, rds, task_id, start_density)
  all_N <- bind_N(all_N, rds, task_id, start_density)
  all_beta_p <- bind_beta_p(all_beta_p, rds, task_id, start_density)
  all_methods <- bind_methods(all_methods, rds, task_id, start_density)
  all_psrf <- bind_psrf(all_psrf, rds, task_id, start_density)

  setTxtProgressBar(pb, i)
}
close(pb)


path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, model_dir, density_dir)
if(!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

write_rds(all_methods, file.path(path, "method_parameter_lookup.rds"))

data_model_summaries <- list(
  # theta = all_theta,
  potential_area = all_area,
  p = all_p
)

write_rds(data_model_summaries, file.path(path, "data_model_summaries.rds"))
message("\ndata model summaries done\n")

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
              med = quantile(value, 0.5),
              high = quantile(value, 0.975))
}


recovery_list <- list()
residual_list <- list()

## capture probability intercepts ------
beta1_long <- all_samples |>
  select_pivot_longer("beta1") |>
  mutate(method_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
         position = 1)

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

recovery_list$beta1 <- recov_beta1(beta1_long, all_psrf)
residual_list$beta1 <- resid_beta1(beta1_long)

message("\ncapture intercepts done\n")

## capture probability covariates ------
beta_p_long <- all_samples |>
  select_pivot_longer("beta_p") |>
  mutate(method_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
         position = as.numeric(str_extract(node, "(?<=\\, )\\d")) + 1)

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

recovery_list$beta_p <- recov_beta_p(beta_p_long, all_psrf)
residual_list$beta_p <- resid_beta_p(beta_p_long)

message("\ncapture covariates done\n")

## gamma ------

gH <- all_methods |>
  select(idx, gamma, method, simulation) |>
  filter(method %in% c("Snares", "Traps")) |>
  mutate(idx = idx - 3) |>
  rename(actual = gamma)

gamma_long <- all_samples |>
  select_pivot_longer("log_gamma[") |>
  mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d"))) |>
  mutate(value = exp(value))

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

recovery_list$gamma <- recov_gamma(gamma_long, gH, all_psrf)
residual_list$gamma <- resid_gamma(gamma_long, gH)

message("\nsaturation constant done\n")

## rho ------
rH <- all_methods |>
  select(idx, rho, method, simulation) |>
  rename(actual = rho)

rho_long<- all_samples |>
  select_pivot_longer("log_rho[") |>
  mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d"))) |>
  mutate(value = exp(value))

recovery_list$rho <- recov_gamma(rho_long, rH, all_psrf)
residual_list$rho <- resid_gamma(rho_long, rH)

message("\nsearch area done\n")

## unique area ------
pH <- all_methods |>
  filter(idx %in% c(1, 4, 5)) |>
  select(idx, p_unique, method, simulation) |>
  rename(actual = p_unique) |>
  mutate(idx = if_else(idx == 1, idx, idx - 2))

p_mu_long <- all_samples |>
  select_pivot_longer("p_mu[") |>
  mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d"))) |>
  mutate(value = ilogit(value))

recovery_list$p_mu <- recov_gamma(p_mu_long, pH, all_psrf)
residual_list$p_mu <- resid_gamma(p_mu_long, pH)

message("\nunique area done\n")

## litter size ------
actual <- 5.290323
ls_long <- all_samples |>
  select_pivot_longer("log_nu") |>
  mutate(value = exp(value))

ls_recovery <- ls_long |>
  group_by(simulation, node, start_density) |>
  my_summary() |>
  mutate(actual = actual) |>
  ungroup() |>
  recovered() |>
  left_join(all_psrf)

ls_residual <- ls_long |>
  mutate(actual = actual) |>
  mutate(value = value - actual) |>
  group_by(node, start_density) |>
  my_summary() |>
  ungroup()

recovery_list$litter_size <- ls_recovery
residual_list$litter_size <- ls_residual

message("\nliter size done\n")

## survival ------
actual <- config$phi_mu
phi_long <- all_samples |>
  select_pivot_longer("phi_mu")

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

recovery_list$phi_mu <- recov_phi(phi_long, actual, all_psrf)
residual_list$phi_mu <- resid_phi(phi_long, actual)

actual <- config$psi_phi
psi_phi_long <- all_samples |>
  select_pivot_longer("psi_phi")

recovery_list$psi_phi <- recov_phi(psi_phi_long, actual, all_psrf)
residual_list$psi_phi <- resid_phi(psi_phi_long, actual)

message("\nsurvival done\n")

write_rds(recovery_list, file.path(path, "parameterRecovery.rds"))
write_rds(residual_list, file.path(path, "parameterResidual.rds"))

## abundance ------
abundance <- all_N |>
  rename(abundance = N)



get_xn <- function(df, H){
  df |>
    select_pivot_longer("N[") |>
    filter(!is.na(value)) |>
    mutate(n_id = as.numeric(str_extract(node, "(?<=\\[)\\d*"))) |>
    left_join(H) |>
    filter(!is.na(abundance)) |>
    mutate(estimated_density = value / property_area)
}

xn <- all_samples |>
  get_xn(abundance)

vals <- c("abundance", "density", "value", "estimated_density")

n_attributes <- xn |>
  select(-value, -estimated_density) |>
  distinct()

abundance_summaries <- xn |>
  select(simulation, property, PPNum, all_of(vals)) |>
  group_by(simulation, property, PPNum) |>
  summarise(low_abundance = quantile(value, 0.025),
            med_abundance = quantile(value, 0.5),
            high_abundance = quantile(value, 0.975),
            var_abundance = var(value),
            low_density = quantile(estimated_density, 0.025),
            med_density = quantile(estimated_density, 0.5),
            high_density = quantile(estimated_density, 0.975),
            var_density = var(estimated_density)) |>
  ungroup() |>
  left_join(n_attributes)

write_rds(abundance_summaries, file.path(path, "abundance_summaries.rds"))
message("\nposterior abundance done\n")

error_by_observation <- xn |>
  select(simulation, property, PPNum, all_of(vals)) |>
  group_by(simulation, property, PPNum) |>
  summarise(mpe_abundance = mean(abs((value+1) - (abundance+1))/(abundance+1))*100,
            mpe_density = mean(abs((estimated_density+0.1) - (density+0.1))/(density+0.1))*100,
            mbias_abundance = mean(value - abundance),
            mbias_density = mean(estimated_density - density),
            mse_abundance = mean((value - abundance)^2),
            mse_density = mean((estimated_density - density)^2),
            rmse_abundance = sqrt(mse_abundance),
            rmse_density = sqrt(mse_density),
            nm_rmse_abundance = rmse_abundance / mean(abundance),
            nm_rmse_density = rmse_density / mean(density)) |>
  ungroup() |>
  arrange(simulation, property, PPNum) |>
  group_by(simulation, property) |>
  mutate(delta = PPNum - lag(PPNum)) |>
  ungroup() |>
  left_join(n_attributes)

write_rds(error_by_observation, file.path(path, "abundance_error_by_observation.rds"))
message("\nabundance error by observation done\n")

error_by_simulation <- xn |>
  select(simulation, all_of(vals)) |>
  group_by(simulation, start_density) |>
  summarise(mpe_abundance = mean(abs((value+1) - (abundance+1))/(abundance+1))*100,
            mpe_density = mean(abs((estimated_density+0.1) - (density+0.1))/(density+0.1))*100,
            mbias_abundance = mean(value - abundance),
            mbias_density = mean(estimated_density - density),
            mse_abundance = mean((value - abundance)^2),
            mse_density = mean((estimated_density - density)^2),
            rmse_abundance = sqrt(mse_abundance),
            rmse_density = sqrt(mse_density),
            nm_rmse_abundance = rmse_abundance / mean(abundance),
            nm_rmse_density = rmse_density / mean(density),
            ns_rmse_abundance = rmse_abundance / sd(abundance),
            ns_rmse_density = rmse_density / sd(density),
            sd_ratio_abundance = sd(value) / sd(abundance),
            sd_ratio_density = sd(estimated_density) / sd(density)) |>
  ungroup()


write_rds(error_by_simulation, file.path(path, "abundance_error_by_simulation.rds"))
message("\nabundance error by simulation done\n")


get_post_take <- function(df, H){
  df |>
    pivot_longer(cols = -c(simulation, start_density),
                 names_to = "p_id") |>
    filter(!is.na(value)) |>
    mutate(p_id = as.numeric(p_id)) |>
    left_join(H)
}

yy <- all_y |>
  get_post_take(all_take)

take_summaries <- yy |>
  group_by(simulation, p_id, start_density) |>
  my_summary() |>
  ungroup() |>
  left_join(all_take)

write_rds(take_summaries, file.path(path, "take_summaries.rds"))
message("\nposterior take error by observation done\n")

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

error_by_observation <- yy |>
  group_by(simulation, p_id, start_density) |>
  take_calc() |>
  ungroup() |>
  select(-nm_rmse, -sd_ratio) |>
  left_join(all_take)

write_rds(error_by_observation, file.path(path, "take_error_by_observation.rds"))
message("\nposterior take error by observation done\n")

error_by_simulation <- yy |>
  group_by(start_density, simulation) |>
  take_calc() |>
  ungroup()

write_rds(error_by_simulation, file.path(path, "take_error_by_simulation.rds"))
message("\nposterior take error by observation done\n")

error_by_simulation_method <- yy |>
  group_by(simulation, method) |>
  take_calc() |>
  ungroup()

write_rds(error_by_simulation_method, file.path(path, "take_error_by_simulation_method.rds"))
message("\nposterior take error by simulation method done\n")
message("=== DONE ===")

