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
density <- a_join |> select(start_density, PPNum, contains("property"), contains("density"), obs_flag) |>
  mutate(recovered = if_else(density >= low_density & density <= high_density, "Recovered", "Not Recovered"),
         extinct = if_else(density == 0, "Extinct", "Extant"))

message("\nProportion recovered by extinct vs extant:")
density |>
  group_by(extinct) |>
  count(recovered) |>
  mutate(proportion = round(n / sum(n), 3))


f_name <- "all_take.rds"
all_take <- map_files2(density_dirs, f_name) |>
  left_join(property_ids) |>
  select(-order, -theta, -p, -property, -obs_flag, -simulation, -p_id)

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

f_name <- "take_effort_summary.rds"
take_effort_summary <- map_files2(density_dirs, f_name) |>
  left_join(property_ids) |>
  ungroup()

density_obs_take <- left_join(take_effort_summary, density_obs)
density_obs_take_params <- left_join(density_obs_take, all_methods, relationship = "many-to-many")

n_methods_prop <- all_take |>
  select(property_id, PPNum, method) |>
  distinct() |>
  count(property_id) |>
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
  left_join(n_methods_prop)

data_final_join <- left_join(density_obs_take_params, n_return)

rescale_variable <- function(x) (x - mean(x)) / sd(x)

data <- data_final_join |>
  filter(density > 0) |>
  pivot_wider(names_from = parameter_name,
              values_from = parameter_value) |>
  mutate(method = as.factor(method),
         methods_used = as.factor(methods_used)) |>
  group_by(method) |>
  mutate(mean_effort = rescale_variable(mean_effort),
         sum_effort = rescale_variable(sum_effort),
         mean_trap_count = rescale_variable(mean_trap_count),
         sum_trap_count = rescale_variable(sum_trap_count),
         n_reps = rescale_variable(n_reps)) |>
  ungroup() |>
  group_by(methods_used) |>
  mutate(return_interval = rescale_variable(return_interval)) |>
  ungroup() |>
  mutate(med_density = rescale_variable(med_density),
         n_methods_used = rescale_variable(n_methods_used),
         sum_take_density = rescale_variable(sum_take_density),
         property_area = rescale_variable(property_area))

path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, "GLMs", model_dir)
if(!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)


## glm to explain normalised rmse given mean effort
p <- 1.01

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

if(task_id == 1){
  message("Fit full model glm")
  mFull <- glmer(nm_rmse_density ~                         # global intercept
                   (1 | method) +                          # random effect intercept for each method
                   (1 | methods_used) +                    # random effect intercept for methods used in each PP
                   (0 + mean_effort | method) +            # mean effort per trap deviations across methods
                   (0 + sum_effort | method) +             # total effort per trap deviations across methods
                   (0 + mean_trap_count | method) +        # mean number of traps deviations across methods
                   (0 + sum_trap_count | method) +         # total number of traps deviations across methods
                   (0 + n_reps | method) +                 # number events deviations across methods
                   (0 + return_interval | methods_used) +  # return interval deviations across levels of methods used
                   med_density +                           # the effect (slope) of median estimated density
                   n_methods_used +                        # the effect (slope) of number of methods used in PP
                   sum_take_density +                      # the effect (slope) of total take in PP as a density
                   property_area +                         # the effect (slope) of property area
                   I(med_density * n_methods_used) +       # the interaction between median estimated density and number of methods used in PP
                   I(med_density * sum_take_density) +     # the interaction between median estimated density and total take in PP as a density
                   I(med_density * property_area) +        # the interaction between median estimated density and property area
                   I(n_methods_used * sum_take_density) +  # the interaction between number of methods used in PP and total take in PP as a density
                   I(n_methods_used * property_area) +     # the interaction between number of methods used in PP property area
                   I(sum_take_density * property_area),    # the interaction between total take in PP as a density and property area
                 family = Tweedie(p = p),
                 data = data)

  summary(mFull)
  write_rds(mFull, file.path(path, "glmer_full.rds"))

} else if(task_id == 2){
  message("Fit full model gam")
  mFullg <- gam(nm_rmse_density ~                         # global intercept
                  s(method, bs = "re") +                          # random effect intercept for each method
                  s(methods_used, bs = "re") +                    # random effect intercept for methods used in each PP
                  s(mean_effort, by = method) +            # mean effort per trap deviations across methods
                  s(sum_effort, by = method) +             # total effort per trap deviations across methods
                  s(mean_trap_count, by = method) +        # mean number of traps deviations across methods
                  s(sum_trap_count, by = method) +         # total number of traps deviations across methods
                  s(n_reps, by = method) +                 # number events deviations across methods
                  s(return_interval, by = methods_used) +  # return interval deviations across levels of methods used
                  s(med_density, bs = "cs") +                           # the effect (slope) of median estimated density
                  s(n_methods_used, bs = "cs") +                        # the effect (slope) of number of methods used in PP
                  s(sum_take_density, bs = "cs") +                      # the effect (slope) of total take in PP as a density
                  s(property_area, bs = "cs") +                         # the effect (slope) of property area
                  s(med_density, n_methods_used) +       # the interaction between median estimated density and number of methods used in PP
                  s(med_density, sum_take_density) +     # the interaction between median estimated density and total take in PP as a density
                  s(med_density, property_area) +        # the interaction between median estimated density and property area
                  s(n_methods_used, sum_take_density) +  # the interaction between number of methods used in PP and total take in PP as a density
                  s(n_methods_used, property_area) +     # the interaction between number of methods used in PP property area
                  s(sum_take_density, property_area),    # the interaction between total take in PP as a density and property area
                family = Tweedie(p = p),
                data = data)

  summary(mFullg)
  write_rds(mFullg, file.path(path, "gam_full.rds"))

}






