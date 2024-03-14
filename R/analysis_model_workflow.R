library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(recipes)
library(rsample)
library(xgboost)

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
  ungroup() |>
  filter(density > 0)

responses <- c("nm_rmse_density", "mpe_density", "mbias_density")

# hyperparameter grid
hyper_grid <- expand_grid(
  responses = responses,
  eta = 0.05,
  max_depth = 3:8,
  min_child_weight = c(0.5, 1, 1.5),
  subsample = 0.5,
  colsample_bytree = c(0.5, 0.75, 1),
  gamma = c(0, 1, 10, 100, 1000),
  lambda = c(0, 1e-2, 0.1, 1, 100, 1000, 10000),
  alpha = c(0, 1e-2, 0.1, 1, 100, 1000, 10000),
  rmse = 0,
  trees = 0
)

n_by_response <- hyper_grid |>
  group_by(responses) |>
  count() |>
  pull(n) |>
  unique()

n_models_per_array <- 210

array_nums_1 <- rep(seq(1, ceiling(n_by_response / n_models_per_array)), each = n_models_per_array)
array_nums_2 <- array_nums_1 + max(array_nums_1)
array_nums_3 <- array_nums_1 + max(array_nums_1)*2
max(array_nums_3)

hyper_grid <- hyper_grid |>
  mutate(array = c(array_nums_1, array_nums_2, array_nums_3))

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
if(is.na(task_id)) task_id <- 6
message("task id: ", task_id)

array_grid <- hyper_grid |>
  filter(array == task_id) |>
  select(-array)

y <- array_grid |>
  pull(responses)

testthat::expect_equal(length(unique(y)), 1)
y <- y[1]

message("\ny: ", y)

df_model <- subset_rename(data, y)
message("\nTraining data:")
glimpse(df_model$train)

baked_data <- my_recipe(df_model$train, df_model$test)

train_data <- baked_data$df_train
X <- train_data |>
  select(-y) |>
  as.matrix()
Y <- train_data |> pull(y)

message("Fitting xgBoost...")
start_time <- Sys.time()

# grid search
for(i in seq_len(nrow(array_grid))) {

  if(i %% 10 == 0){
    message("[", i, "/", nrow(array_grid), "] ", round(i/nrow(array_grid)*100), "%")
  }

  set.seed(123)
  m <- xgb.cv(
    data = X,
    label = Y,
    nrounds = 5000,
    objective = "reg:squarederror",
    metrics = "rmse",
    early_stopping_rounds = 50,
    nfold = 10,
    verbose = 0,
    params = list(
      eta = array_grid$eta[i],
      max_depth = array_grid$max_depth[i],
      min_child_weight = array_grid$min_child_weight[i],
      subsample = array_grid$subsample[i],
      colsample_bytree = array_grid$colsample_bytree[i],
      gamma = array_grid$gamma[i],
      lambda = array_grid$lambda[i],
      alpha = array_grid$alpha[i]
    )
  )
  array_grid$rmse[i] <- min(m$evaluation_log$test_rmse_mean)
  array_grid$trees[i] <- m$best_iteration
}
message("xgBoost complete!")

out <- array_grid |>
  as_tibble() |>
  arrange(rmse)

message("All fits")
print(out)

path <- file.path(top_dir, project_dir, analysis_dir, dev_dir, "gradientBoosting")
filename <- file.path(path, paste0("xgbTree_", task_id, ".rds"))
write_rds(out, filename)

total_time <- Sys.time() - start_time
message("Elapsed time: ")
print(total_time)

message("=== DONE ===")


