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

take_joint_return <- left_join(all_take, n_return)

take_effort_summary <- take_joint_return |>
  group_by(property_id, PPNum, methods_used) |>
  summarise(mean_effort_per_unit = mean(effort_per),
            sum_effort_per_unit = sum(effort_per),
            mean_effort = mean(effort),
            sum_effort = sum(effort),
            mean_unit_count = mean(trap_count),
            sum_unit_count = sum(trap_count),
            n_reps_pp = n()) |>
  ungroup()

data_final_join <- left_join(density_obs, take_effort_summary)

rescale_variable <- function(x) (x - mean(x)) / (2 * sd(x))

outlier <- data_final_join |>
  filter(density > 0) |>
  pull(nm_rmse_abundance) |>
  quantile(0.995)

data <- data_final_join |>
  ungroup() |>
  filter(density > 0,
         nm_rmse_abundance < outlier) |>
  select(property_id, PPNum, property_area, med_density,
         nm_rmse_density, mbias_density, mpe_density,
         sum_take_density, delta, methods_used, mean_effort_per_unit,
         sum_effort_per_unit, mean_effort, sum_effort, mean_unit_count,
         sum_unit_count, n_reps_pp) |>
  distinct() |>
  mutate(property_area = rescale_variable(property_area),
         med_density = rescale_variable(med_density),
         total_take_density = rescale_variable(sum_take_density),
         delta = rescale_variable(delta),
         methods_used = as.factor(methods_used),
         mean_effort_per_unit = rescale_variable(mean_effort_per_unit),
         sum_effort_per_unit = rescale_variable(sum_effort_per_unit),
         mean_effort_raw = rescale_variable(mean_effort),
         sum_effort_raw = rescale_variable(sum_effort),
         mean_unit_count = rescale_variable(mean_unit_count),
         sum_unit_count = rescale_variable(sum_unit_count),
         n_reps_pp = rescale_variable(n_reps_pp)) |>
  select(-mean_effort, -sum_effort, -sum_take_density)




path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, "GLMs", model_dir)
if(!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

models <- expand_grid(
  y = c("nrmse", "bias", "mpe"),
  effort = c("per_unit", "raw"),
  agg = c("mean", "sum")
)

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

source("R/functions_analysis.R")

model_to_run <- models |> slice(task_id)

y <- pull(model_to_run, y)
effort <- pull(model_to_run, effort)
agg <- pull(model_to_run, agg)

fit <- fit_glm_all(
  df = data,
  y = y,
  effort = effort,
  agg = agg,
  path = path
)
warnings()
fit

#  prevent fitting sub-models to different datasets
oop <- options(na.action = "na.fail")

dd <- MuMIn::dredge(fit)
filename <- paste(y, effort, agg, sep = "-")
outname <- file.path(path, paste0(filename, "-dredge.rds"))
write_rds(fit, outname)

dd
