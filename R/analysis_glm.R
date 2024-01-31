library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(mgcv)
library(lme4)

# analysis_dir <- "analysis"
# model_dir <- "betaSurvival_uniqueAreaTrapSnare"
# path <- file.path(analysis_dir, model_dir)
# density_dirs <- list.files(path)

config_name <- "hpc_production"
config <- config::get(config = config_name)
top_dir <- config$top_dir
out_dir <- config$out_dir
analysis_dir <- config$analysis_dir
dev_dir <- config$dev_dir
model_dir <- config$model_dir
project_dir <- config$project_dir
path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, model_dir)
density_dirs <- list.files(path)

density_dirs <- grep("density_", density_dirs, value = TRUE)

map_files <- function(dirs_vec, file_name, node){

  get_files <- function(density_dir, file_name, node){
    sim_results <- file.path(path, density_dir)
    ls <- read_rds(file.path(sim_results, file_name))
    ls[[node]]
  }

  dirs_vec |>
    map(\(x) get_files(x, file_name, node)) |>
    list_rbind() |>
    mutate(start_density = as.factor(start_density))

}

map_files2 <- function(dirs_vec, file_name){

  get_files <- function(density_dir, file_name, node){
    sim_results <- file.path(path, density_dir)
    ls <- read_rds(file.path(sim_results, file_name))
  }

  dirs_vec |>
    map(\(x) get_files(x, file_name)) |>
    list_rbind() |>
    mutate(start_density = as.factor(start_density))

}

message("Get abundance summaries")
abundance_file <- "abundance_summaries.rds"
f_name <- abundance_file
df <- map_files2(density_dirs, f_name)

# create property ID for easier joining
property_ids <- df |>
  select(start_density, simulation, property) |>
  distinct() |>
  mutate(property_id = paste(start_density, simulation, property, sep = "-"))

message("\nSample size by start density:")
sample_size <- df |>
  select(start_density, simulation) |>
  distinct() |>
  group_by(start_density) |>
  count() |>
  rename(`Start density` = start_density,
         `n simulations` = n)
print(sample_size)

abundance_file <- "abundance_error_by_observation.rds"
f_name <- abundance_file
abundance_error_by_observation <- map_files2(density_dirs, f_name) |>
  left_join(property_ids) |>
  mutate(nm_rmse_abundance = if_else(rmse_abundance == 0, 0, nm_rmse_abundance),
         nm_rmse_density = if_else(rmse_density == 0, 0, nm_rmse_density))

abundance_summary <- map_files2(density_dirs, "abundance_summaries.rds") |>
  left_join(property_ids)

a_join <- left_join(abundance_summary, abundance_error_by_observation)
density <- a_join |> select(start_density, PPNum, contains("property"), contains("abundance"), contains("density"), obs_flag) |>
  mutate(recovered = if_else(abundance >= low_abundance & density <= high_abundance, "Recovered", "Not Recovered"),
         extinct = if_else(abundance == 0, "Extinct", "Extant"))

message("\nProportion recovered by extinct vs extant:")
density |>
  group_by(extinct) |>
  count(recovered) |>
  mutate(proportion = round(n / sum(n), 3))


f_name <- "all_take.rds"
all_take <- map_files2(density_dirs, f_name) |>
  left_join(property_ids) |>
  select(-theta, -p, -property, -obs_flag, -simulation, -p_id) |>
  mutate(effort = effort_per * trap_count)

sum_take <- all_take |>
  select(property_id, PPNum, sum_take) |>
  distinct()

density_obs <- density |>
  filter(obs_flag == 1) |>
  left_join(sum_take) |>
  mutate(sum_take_density = sum_take / property_area) |>
  group_by(property_id) |>
  mutate(delta = c(0, diff(PPNum)))

f_name <- "method_parameter_lookup.rds"
all_methods <- map_files2(density_dirs, f_name)  |>
  pivot_longer(cols = c(p_unique, rho, gamma),
               names_to = "parameter_name",
               values_to = "parameter_value") |>
  filter(!is.na(parameter_value))

take_effort_summary <- all_take |>
  group_by(property_id, PPNum, method) |>
  summarise(mean_effort_per = mean(effort_per),
            sum_effort_per = sum(effort_per),
            mean_effort = mean(effort),
            sum_effort = sum(effort),
            mean_trap_count = mean(trap_count),
            sum_trap_count = sum(trap_count),
            n_reps = n())

density_obs_take <- left_join(take_effort_summary, density_obs)
density_obs_take_params <- left_join(density_obs_take, all_methods, relationship = "many-to-many")

n_methods_pp <- all_take |>
  select(property_id, PPNum, method) |>
  distinct() |>
  group_by(property_id, PPNum) |>
  count() |>
  rename(n_methods_used = n)

n_return <- all_take |>
  select(property_id, PPNum, method) |>
  distinct() |>
  pivot_wider(names_from = method,
              values_from = method) |>
  unite(method,
        -c(property_id, PPNum),
        sep = ", ",
        na.rm = TRUE) |>
  group_by(property_id, method) |>
  mutate(return_interval = c(0, diff(PPNum))) |>
  ungroup() |>
  rename(methods_used = method) |>
  left_join(n_methods_pp)

data_final_join <- left_join(density_obs_take, n_return)

rescale_variable <- function(x) (x - mean(x)) / sd(x)

data <- data_final_join |>
  ungroup() |>
  filter(density > 0,
         nm_rmse_abundance < quantile(data$nm_rmse_abundance, 0.999)) |>
  rename(nrmse = nm_rmse_abundance) |>
  select(-abundance, -nm_rmse_density, -rmse_abundance, -rmse_density) |>
  mutate(method = as.factor(method),
         methods_used = as.factor(methods_used),
         property_id = as.factor(property_id)) |>
  group_by(method) |>
  mutate(mean_effort_method = rescale_variable(mean_effort),
         sum_effort_method = rescale_variable(sum_effort),
         mean_effort_per_method = rescale_variable(mean_effort_per),
         sum_effort_per_method = rescale_variable(sum_effort_per),
         mean_trap_count_method = rescale_variable(mean_trap_count),
         sum_trap_count_method = rescale_variable(sum_trap_count),
         n_reps_method = rescale_variable(n_reps)) |>
  ungroup() |>
  # group_by(methods_used) |>
  # mutate() |>
  # ungroup() |>
  mutate(med_density = rescale_variable(med_density),
         return_interval = rescale_variable(return_interval),
         #n_methods_used = rescale_variable(n_methods_used),
         sum_take_density = rescale_variable(sum_take_density),
         property_area = rescale_variable(property_area),
         start_density = as.factor(start_density),
         mean_effort = rescale_variable(mean_effort),
         sum_effort = rescale_variable(sum_effort),
         mean_effort_per = rescale_variable(mean_effort_per),
         sum_effort_per = rescale_variable(sum_effort_per),
         mean_trap_count = rescale_variable(mean_trap_count),
         sum_trap_count = rescale_variable(sum_trap_count),
         n_reps = rescale_variable(n_reps))

path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, "GLMs", model_dir)
if(!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)


args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])


if(task_id == 0){

  ## Null model - random intercepts by methods used in PP
  message("\n=== Model 0 ===")
  message("  GLM")
  m_list <- list()
  glm <- glmer(nrmse ~
                   (1 | methods_used),
                 family = Gamma(link = "log"),
                 data = data)
  m_list$glm <- glm
  write_rds(m_list, file.path(path, paste0("m", task_id, ".rds")))
  message("    -> done")

  message("  GAM")
  gam <- gam(nrmse ~
                 s(methods_used, bs = "re"),
               family = Gamma(link = "log"),
               data = data)
  m_list$m0 <- gam
  write_rds(m_list, file.path(path, paste0("m", task_id, ".rds")))
  message("    -> done")

} else if(task_id == 1){

  ## All individual effects +
  ## mean effort and mean trap count in a PP across methods
  message("\n=== Model 1 ===")
  message("  GLM")
  m_list <- list()
  glm <- glmer(nrmse ~
                   (1 | methods_used) +
                   return_interval +
                   med_density +
                   delta +
                   n_reps +
                   n_methods_used +
                   sum_take_density +
                   property_area +
                   mean_effort +
                   mean_trap_count,
                 family = Gamma(link = "log"),
                 data = data)
  m_list$glm <- glm
  write_rds(m_list, file.path(path, paste0("m", task_id, ".rds")))
  message("    -> done")

  message("  GAM")
  gam <- gam(nrmse ~
                 s(methods_used, bs = "re") +
                 s(return_interval, bs = "cr") +
                 s(med_density, bs = "cr") +
                 s(delta, bs = "cr") +
                 s(n_reps, bs = "cr") +
                 s(n_methods_used, bs = "cr") +
                 s(sum_take_density, bs = "cr") +
                 s(property_area, bs = "cr") +
                 s(mean_effort, bs = "cr") +
                 s(mean_trap_count, bs = "cr"),
               family = Gamma(link = "log"),
               data = data)
  m_list$glm <- gam
  write_rds(m_list, file.path(path, paste0("m", task_id, ".rds")))
  message("    -> done")

} else if(task_id == 2){

  ## All individual effects +
  ## sum effort and sum trap count in a PP across methods
  message("\n=== Model 2 ===")
  message("  GLM")
  m_list <- list()
  glm <- glmer(nrmse ~
                 (1 | methods_used) +
                 return_interval +
                 med_density +
                 delta +
                 n_reps +
                 n_methods_used +
                 sum_take_density +
                 property_area +
                 sum_effort +
                 sum_trap_count,
               family = Gamma(link = "log"),
               data = data)
  m_list$glm <- glm
  write_rds(m_list, file.path(path, paste0("m", task_id, ".rds")))
  message("    -> done")

  message("  GAM")
  gam <- gam(nrmse ~
               s(methods_used, bs = "re") +
               s(return_interval, bs = "cr") +
               s(med_density, bs = "cr") +
               s(delta, bs = "cr") +
               s(n_reps, bs = "cr") +
               s(n_methods_used, bs = "cr") +
               s(sum_take_density, bs = "cr") +
               s(property_area, bs = "cr") +
               s(sum_effort, bs = "cr") +
               s(sum_trap_count, bs = "cr"),
             family = Gamma(link = "log"),
             data = data)
  m_list$glm <- gam
  write_rds(m_list, file.path(path, paste0("m", task_id, ".rds")))
  message("    -> done")

} else if(task_id == 3){

  ## All individual effects +
  ## mean effort and mean trap count in a PP for each method +
  ## random intercept by method
  message("\n=== Model 3 ===")
  message("  GLM")
  m_list <- list()
  glm <- glmer(nrmse ~
                 (1 | methods_used) +
                 (1 | method) +
                 return_interval +
                 med_density +
                 delta +
                 n_reps +
                 n_methods_used +
                 sum_take_density +
                 property_area +
                 mean_effort_method +
                 mean_trap_count_method,
               family = Gamma(link = "log"),
               data = data)
  m_list$glm <- glm
  write_rds(m_list, file.path(path, paste0("m", task_id, ".rds")))
  message("    -> done")

  message("  GAM")
  gam <- gam(nrmse ~
               s(methods_used, bs = "re") +
               s(method, bs = "re") +
               s(return_interval, bs = "cr") +
               s(med_density, bs = "cr") +
               s(delta, bs = "cr") +
               s(n_reps, bs = "cr") +
               s(n_methods_used, bs = "cr") +
               s(sum_take_density, bs = "cr") +
               s(property_area, bs = "cr") +
               s(mean_effort_method, bs = "cr") +
               s(mean_trap_count_method, bs = "cr"),
             family = Gamma(link = "log"),
             data = data)
  m_list$glm <- gam
  write_rds(m_list, file.path(path, paste0("m", task_id, ".rds")))
  message("    -> done")

} else if(task_id == 4){

  ## All individual effects +
  ## sum effort and sum trap count in a PP for each method +
  ## random intercept by method
  message("\n=== Model 4 ===")
  message("  GLM")
  m_list <- list()
  glm <- glmer(nrmse ~
                 (1 | methods_used) +
                 (1 | method) +
                 return_interval +
                 med_density +
                 delta +
                 n_reps +
                 n_methods_used +
                 sum_take_density +
                 property_area +
                 sum_effort_method +
                 sum_trap_count_method,
               family = Gamma(link = "log"),
               data = data)
  m_list$glm <- glm
  write_rds(m_list, file.path(path, paste0("m", task_id, ".rds")))
  message("    -> done")

  message("  GAM")
  gam <- gam(nrmse ~
               s(methods_used, bs = "re") +
               s(method, bs = "re") +
               s(return_interval, bs = "cr") +
               s(med_density, bs = "cr") +
               s(delta, bs = "cr") +
               s(n_reps, bs = "cr") +
               s(n_methods_used, bs = "cr") +
               s(sum_take_density, bs = "cr") +
               s(property_area, bs = "cr") +
               s(sum_effort_method, bs = "cr") +
               s(sum_trap_count_method, bs = "cr"),
             family = Gamma(link = "log"),
             data = data)
  m_list$glm <- gam
  write_rds(m_list, file.path(path, paste0("m", task_id, ".rds")))
  message("    -> done")

}


