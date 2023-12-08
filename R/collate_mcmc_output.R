# --------------------------------------------------------------------
#
# Script for collating simulation results
#
# John Foster
#
# --------------------------------------------------------------------

library(tidyverse)
library(nimble)
library(coda)

config <- config::get(config = "default")

out_dir <- config$out_dir
model_dir <- config$model_dir

sim_results <- file.path(out_dir, model_dir)
tasks <- 11

add_ids <- function(df, task_id, s_density){
  df |>
    mutate(simulation = task_id,
           start_density = s_density)
}

all_samples <- tibble()
all_take <- tibble()
all_N <- tibble()
all_beta_p <- tibble()
all_methods <- tibble()
all_y <- tibble()
all_area <- tibble()

for(i in tasks){
  task_dir <- file.path(sim_results, i)
  rds <- read_rds(file.path(task_dir, "simulation_data.rds"))

  bad_mcmc <- rds$bad_mcmc #| any(rds$psrf > 1.3)
  task_id <- i
  start_density <- rds$start_density
  # already_collated <- task_id %in% prev_tasks

  # if(bad_mcmc | already_collated) next

  samples <- rds$posterior_samples |>
    add_ids(task_id, start_density)

  all_samples <- bind_rows(all_samples, samples)

  y_pred <- rds$posterior_take
  colnames(y_pred) <- 1:ncol(y_pred)
  y_pred <- y_pred |>
    as_tibble() |>
    add_ids(task_id, start_density)
  all_y <- bind_rows(all_y, y_pred)

  pot_area <- rds$posterior_potential_area
  pot_area <- as_tibble(pot_area) |>
    add_ids(task_id, start_density) |>
    mutate(p_id = 1:n())
  all_area <- bind_rows(all_area, pot_area)

  take <- rds$take |>
    add_ids(task_id, start_density) |>
    mutate(p_id = 1:n())
  all_take <- bind_rows(all_take, take)

  obs_flag <- take |>
    select(property, county, PPNum) |>
    distinct() |>
    mutate(obs_flag = 1)

  N <- rds$N |>
    add_ids(task_id, start_density) |>
    left_join(obs_flag) |>
    mutate(obs_flag = if_else(is.na(obs_flag), 0, obs_flag))
  all_N <- bind_rows(all_N, N)

  # need a lookup table for known data model covariates
  bH <- tibble(
    method_idx = rep(1:nrow(rds$beta_p), ncol(rds$beta_p)),
    position = rep(1:ncol(rds$beta_p), each = nrow(rds$beta_p)),
    actual = as.numeric(rds$beta_p)
  ) |>
    add_ids(task_id, start_density)

  all_beta_p <- bind_rows(all_beta_p, bH)

  method_lookup <- rds$method_lookup |>
    add_ids(task_id, start_density)

  all_methods <- bind_rows(all_methods, method_lookup)

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
              med = quantile(value, 0.5),
              high = quantile(value, 0.975))
}


path <- file.path("analysis", model_dir)
if(!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

recovery_list <- list()
residual_list <- list()

## capture probability intercepts ------
beta1_long <- all_samples |>
  select_pivot_longer("beta1") |>
  mutate(method_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
         position = 1)

beta1_recovery <- beta1_long |>
  group_by(simulation, node, method_idx, position, start_density) |>
  my_summary() |>
  left_join(all_beta_p) |>
  ungroup() |>
  recovered()

beta1_residual <- beta1_long |>
  left_join(all_beta_p) |>
  mutate(value = value - actual) |>
  group_by(node, position, method_idx, start_density) |>
  my_summary() |>
  ungroup()

recovery_list$beta1 <- beta1_recovery
residual_list$beta1 <- beta1_residual

message("capture intercepts done")

## capture probability covariates ------
beta_p_long <- all_samples |>
  select_pivot_longer("beta_p") |>
  mutate(method_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
         position = as.numeric(str_extract(node, "(?<=\\, )\\d")) + 1)

beta_p_recovery <- beta_p_long |>
  group_by(simulation, node, method_idx, position, start_density) |>
  my_summary() |>
  left_join(all_beta_p) |>
  ungroup() |>
  recovered()

beta_p_residual <- beta_p_long |>
  left_join(all_beta_p) |>
  mutate(value = value - actual) |>
  group_by(node, position, method_idx, start_density) |>
  my_summary() |>
  ungroup()

recovery_list$beta_p <- beta_p_recovery
residual_list$beta_p <- beta_p_residual

message("capture covariates done")

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

gamma_recovery <- gamma_long |>
  group_by(simulation, node, idx, start_density) |>
  my_summary() |>
  left_join(gH) |>
  ungroup() |>
  recovered()

gamma_residual <- gamma_long |>
  left_join(gH)|>
  mutate(value = value - actual) |>
  group_by(node, idx, start_density) |>
  my_summary() |>
  ungroup()

recovery_list$gamma <- gamma_recovery
residual_list$gamma <- gamma_residual

message("saturation constant done")

## rho ------
rH <- all_methods |>
  select(idx, rho, method, simulation) |>
  rename(actual = rho)

rho_long<- all_samples |>
  select_pivot_longer("log_rho[") |>
  mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d"))) |>
  mutate(value = exp(value))

rho_recovery <- rho_long |>
  group_by(simulation, node, idx, start_density) |>
  my_summary() |>
  left_join(rH) |>
  ungroup() |>
  recovered()

rho_residual <- rho_long |>
  left_join(rH)|>
  mutate(value = value - actual) |>
  group_by(node, idx, start_density) |>
  my_summary() |>
  ungroup()

recovery_list$rho <- rho_recovery
residual_list$rho <- rho_residual

message("search area done")

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

p_mu_recovery <- p_mu_long |>
  group_by(simulation, node, idx, start_density) |>
  my_summary() |>
  left_join(pH) |>
  ungroup() |>
  recovered()

p_mu_residual <- p_mu_long |>
  left_join(pH)|>
  mutate(value = value - actual) |>
  group_by(node, idx, start_density) |>
  my_summary() |>
  ungroup()

recovery_list$p_mu <- p_mu_recovery
residual_list$p_mu <- p_mu_residual

message("unique area done")

## litter size ------
actual <- 5.290323
ls_long <- all_samples |>
  select_pivot_longer("log_mean_ls") |>
  mutate(value = exp(value))

ls_recovery <- ls_long |>
  group_by(simulation, node, start_density) |>
  my_summary() |>
  mutate(actual = actual) |>
  ungroup() |>
  recovered()

ls_residual <- ls_long |>
  mutate(actual = actual) |>
  mutate(value = value - actual) |>
  group_by(node, start_density) |>
  my_summary() |>
  ungroup()

recovery_list$litter_size <- ls_recovery
residual_list$litter_size <- ls_residual

message("liter size done")

## survival ------
actual <- config$phi_mu
phi_long <- all_samples |>
  select_pivot_longer("phi_mu")

phi_recovery <- phi_long |>
  group_by(simulation, node, start_density) |>
  my_summary() |>
  mutate(actual = actual) |>
  ungroup() |>
  recovered()

phi_residual <- phi_long |>
  mutate(actual = actual) |>
  mutate(value = value - actual) |>
  group_by(node, start_density) |>
  my_summary() |>
  ungroup()

recovery_list$phi_mu <- phi_recovery
residual_list$phi_mu <- phi_residual

actual <- config$psi_phi
psi_phi_long <- all_samples |>
  select_pivot_longer("psi_phi")

psi_phi_recovery <- psi_phi_long |>
  group_by(simulation, node, start_density) |>
  my_summary() |>
  mutate(actual = actual) |>
  ungroup() |>
  recovered()

psi_phi_residual <- psi_phi_long |>
  mutate(actual = actual) |>
  mutate(value = value - actual) |>
  group_by(node, start_density) |>
  my_summary() |>
  ungroup()

recovery_list$psi_phi <- psi_phi_recovery
residual_list$psi_phi <- psi_phi_residual

message("survival done")
write_rds(recovery_list, file.path(path, "parameterRecovery.rds"))
write_rds(residual_list, file.path(path, "parameterResidual.rds"))

## abundance ------
abundance <- all_N |>
  rename(abundance = N)

xn <- all_samples |>
  select_pivot_longer("xn[") |>
  filter(!is.na(value)) |>
  mutate(n_id = as.numeric(str_extract(node, "(?<=\\[)\\d*"))) |>
  left_join(abundance) |>
  filter(!is.na(abundance)) |>
  mutate(estimated_density = value / property_area)

xn_posterior <- xn |>
  group_by(node, n_id, county, property, PPNum, abundance, simulation, property_area, density, start_density, obs_flag) |>
  summarise(low_abundance = quantile(value, 0.025),
            med_abundance = quantile(value, 0.5),
            high_abundance = quantile(value, 0.975),
            var_abundnace = var(value),
            low_density = quantile(estimated_density, 0.025),
            med_density = quantile(estimated_density, 0.5),
            high_density = quantile(estimated_density, 0.975),
            var_density = var(estimated_density)) |>
  ungroup()

message("posterior abundance done")

xn_error <- xn |>
  group_by(node, n_id, county, property, PPNum, abundance, simulation, property_area, density, start_density, obs_flag) |>
  summarise(mae_abundance = mean(abs(value - abundance)),
            mae_density = mean(abs(estimated_density - density)),
            mpe_abundance = mean(abs((value+1) - (abundance+1))/(abundance+1))*100,
            mpe_density = mean(abs((estimated_density+0.1) - (density+0.1))/(density+0.1))*100,
            mbias_abundance = mean(value - abundance),
            mbias_density = mean(estimated_density - density),
            mse_abundance = mean((value - abundance)^2),
            mse_density = mean((estimated_density - density)^2),
            rmse_abundance = (sqrt(mse_abundance)),
            rmse_density = (sqrt(mse_density))) |>
  ungroup() |>
  arrange(simulation, property, PPNum) |>
  group_by(simulation, property) |>
  mutate(delta = PPNum - lag(PPNum)) |>
  ungroup()


abundance_list <- list(
  abundance_summaries = xn_posterior,
  abundance_metrics = xn_error
)
write_rds(abundance_list, file.path(path, "abundance.rds"))
message("abundance metrics done")

posterior_take <- all_y |>
  pivot_longer(cols = -c(simulation, start_density),
               names_to = "p_id") |>
  filter(!is.na(value)) |>
  mutate(p_id = as.numeric(p_id)) |>
  left_join(all_take)

take_summaries <- posterior_take |>
  group_by(simulation, p_id, start_density) |>
  my_summary() |>
  ungroup() |>
  left_join(all_take)

take_metrics <- posterior_take |>
  group_by(simulation, p_id, start_density) |>
  summarise(mae_abundance = mean(abs(value - take)),
            mpe_abundance = mean(abs((value+1) - (take+1))/(take+1))*100,
            mbias_abundance = mean(value - take),
            mse_abundance = mean((value - take)^2),
            rmse_abundance = (sqrt(mse_abundance))) |>
  ungroup() |>
  left_join(all_take)

take_list <- list(
  take_summaries = take_summaries,
  take_metrics = take_metrics
)
write_rds(take_list, file.path(path, "take.rds"))
message("posterior take done")

