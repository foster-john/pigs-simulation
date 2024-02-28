library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(recipes)
library(rsample)
library(caret)
library(ranger)

source("R/functions_analysis.R")

analysis_dir <- "analysis"
model_dir <- "betaSurvival_uniqueAreaTrapSnare"
path <- file.path(analysis_dir, model_dir)

config_name <- "hpc_production"
config <- config::get(config = config_name)
top_dir <- config$top_dir
analysis_dir <- config$analysis_dir
dev_dir <- config$dev_dir
model_dir <- config$model_dir
project_dir <- config$project_dir
path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, model_dir)

data <- read_rds(file.path(path, "abundanceScoresByPrimaryPeriod.rds")) |>
  filter(density > 0)
#glimpse(data)

tasks <- expand_grid(
 y = c("rmsle_density", "nm_rmse_density", "mpe_density", "mbias_density"),
 ml = c("ranger", "knn")
)

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
#if(is.na(task_id)) task_id <- 6
message("task id: ", task_id)

y <- tasks |> slice(task_id) |> pull(y)
ml <- tasks |> slice(task_id) |> pull(ml)

message("\ny: ", y)
message("ML: ", ml)

df <- subset_rename(data, y)

samps <- sample.int(nrow(df), 5000, replace = FALSE)
df_model <- df |> slice(samps)
glimpse(df_model)

start_time <- Sys.time()

filename <- file.path(path, paste0(y, "_", ml, ".rds"))
fit_ml(df_model, ml, filename)

total_time <- Sys.time() - start_time
message("Elapsed time: ")
print(total_time)
message("=== DONE ===")

