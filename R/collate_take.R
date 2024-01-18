# --------------------------------------------------------------------
#
# Script for collating take from simulation results
#
# John Foster
#
# --------------------------------------------------------------------

message("=== TAKE ===")

# will get config_name from here
source("R/functions_collate.R")
config <- config::get(config = config_name)

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])

read_path <- get_path("read", config_name, task_id)

density_tasks <- list.files(read_path)
message("Tasks to collate ", length(density_tasks))

nodes_vec <- c(
  "all_y",
  "all_take"
)
tasks_ls <- get_tasks(density_tasks, read_path, nodes_vec)
all_y <- tasks_ls$all_y
all_take <- tasks_ls$all_take

path <- get_path("write", config_name, task_id)
if(!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

yy <- all_y |>
  get_post_take(all_take)

take_summaries <- yy |>
  group_by(simulation, p_id, start_density) |>
  my_summary() |>
  ungroup() |>
  left_join(all_take)

write_rds(take_summaries, file.path(path, "take_summaries.rds"))
rm(take_summaries)
gc()
message("\nposterior take error by observation done\n")

error_by_observation <- yy |>
  group_by(simulation, p_id, start_density) |>
  take_calc() |>
  ungroup() |>
  select(-nm_rmse, -sd_ratio) |>
  left_join(all_take)

write_rds(error_by_observation, file.path(path, "take_error_by_observation.rds"))
rm(error_by_observation)
gc()
message("\nposterior take error by observation done\n")

error_by_simulation <- yy |>
  group_by(start_density, simulation) |>
  take_calc() |>
  ungroup()

write_rds(error_by_simulation, file.path(path, "take_error_by_simulation.rds"))
rm(error_by_simulation)
gc()
message("\nposterior take error by observation done\n")

error_by_simulation_method <- yy |>
  group_by(simulation, start_density, method) |>
  take_calc() |>
  ungroup()

write_rds(error_by_simulation_method, file.path(path, "take_error_by_simulation_method.rds"))
rm(error_by_simulation_method)
gc()
message("\nposterior take error by simulation method done\n")

error_by_property <- yy |>
  group_by(simulation, start_density, property) |>
  take_calc() |>
  ungroup()

write_rds(error_by_property, file.path(path, "take_error_by_property.rds"))
rm(error_by_property)
gc()
message("\nposterior take error by property done\n")



message("=== TAKE DONE ===")

