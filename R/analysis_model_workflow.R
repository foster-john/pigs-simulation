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
  filter(density > 0) |>
  ungroup() |>
  group_by()

# samps <- sample.int(nrow(data), 5000, replace = FALSE)
# data <- data |> slice(samps)

responses <- c("nm_rmse_density", "mpe_density", "mbias_density")

# hyperparameter grid
hyper_grid <- expand_grid(
  responses = c("nm_rmse_density", "mpe_density", "mbias_density"),
  nrounds = c(50, 100, 500, 1000, 2000),
  eta = c(0.05, 0.1),
  lambda = c(0, 1e-2, 0.1, 1, 100, 1000, 10000),
  alpha = c(0, 1e-2, 0.1, 1, 100, 1000, 10000)
)

n_by_response <- hyper_grid |>
  group_by(responses) |>
  count() |>
  pull(n) |>
  unique()

n_models_per_array <- 10

array_nums <- rep(seq(1, ceiling(n_by_response / n_models_per_array)), each = n_models_per_array)
array_nums_2 <- array_nums + max(array_nums)
array_nums_3 <- array_nums + max(array_nums)*2

hyper_grid <- hyper_grid |>
  mutate(array = c(array_nums, array_nums_2, array_nums_3))

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
if(is.na(task_id)) task_id <- 6
message("task id: ", task_id)

y <- hyper_grid |>
  filter(array == task_id) |>
  pull(responses)

testthat::expect_equal(length(unique(y)), 1)
y <- y[1]

message("\ny: ", y)

df_model <- subset_rename(data, y)
glimpse(df_model)

train <- df_model$train
test <- df_model$test

blueprint <- recipe(y ~ ., data = train) |>
  step_dummy(all_nominal_predictors()) |>
  step_nzv(all_predictors()) |>
  step_center(all_numeric_predictors()) |>
  step_scale(all_numeric_predictors())

cv <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 5
)

tune_grid <- hyper_grid |>
  filter(array == task_id) |>
  select(-responses, -array) |>
  as.data.frame()

message("Fitting xgBoost...")
start_time <- Sys.time()

message("  make cluster")
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

fit <- train(
  blueprint,
  data = train,
  method = "xgbLinear",
  trControl = cv,
  tuneGrid = tune_grid,
  metric = "RMSE"
)

stopCluster(cl)
message("  complete!")

total_time <- Sys.time() - start_time
message("Elapsed time: ")
print(total_time)

path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, "gradientBoosting")
filename <- file.path(path, paste0(task_id, ".rds"))

out <- as_tibble(fit$results)
write_rds(out, filename)

message("All fits")
print(out)

message("Best tune")
print(fit$bestTune)

message("=== DONE ===")


