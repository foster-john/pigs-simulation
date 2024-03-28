library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(recipes)
library(rsample)
library(xgboost)
library(doParallel)
library(foreach)
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
  filter(density > 0) |>
  mutate(mbias_density_class = as.numeric(mbias_density > 0)) |>
  rename(mbias_density_reg = mbias_density)

responses <- c("nm_rmse_density", "mpe_density", "mbias_density_reg", "mbias_density_class")

eta_grid <- tibble(
  responses = responses,
  eta = c(0.1, 0.1, 0.3, 0.1))

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
if(is.na(task_id)) task_id <- 1
message("task id: ", task_id)

task <- eta_grid |>
  dplyr::slice(task_id)

y <- task$responses
message("\ny: ", y)

eta <- task$eta
message("eta: ", eta)
message("Number of threads: ", Sys.getenv("OMP_NUM_THREADS"))

df_model <- subset_rename(data, y)
baked_data <- my_recipe(df_model$train, df_model$test)

train_data <- baked_data$df_train
X <- train_data |>
  select(-y) |>
  as.matrix()
Y <- train_data |> pull(y)

objective <- if_else(y == "mbias_density_class", "binary:logistic", "reg:squarederror")

fit_xgBoost <- function(eta, X, Y, objective, omp_thread = Sys.getenv("OMP_NUM_THREADS")){
  m <- xgb.cv(
    data = X,
    label = Y,
    nrounds = 5000,
    objective = objective,
    metrics = "rmse",
    early_stopping_rounds = 50,
    nfold = 10,
    nthread = as.numeric(omp_thread),
    verbose = 0,
    params = list(
      eta = eta
    )
  )
}

timing <- system.time(
  fit_xgBoost(eta, X, Y, objective)
)[3]

timing


stop()


# hyperparameter grid
hyper_grid <- expand_grid(
  responses = responses,
  max_depth = 3:8,
  min_child_weight = c(0.5, 1),
  subsample = c(0.5, 1),
  colsample_bytree = c(0.5, 0.75, 1),
  gamma = c(0, 1, 10, 100, 1000),
  lambda = c(0, 1e-2, 0.1, 1, 100, 1000),
  alpha = c(0, 1e-2, 0.1, 1, 100, 1000),
  rmse = 0,
  trees = 0
)

# grid search
start_time <- Sys.time()
# for(i in seq_len(nrow(array_grid))) {
for(i in 11:20) {

  if(i %% 2 == 0){
    message("[", i, "/", nrow(array_grid), "] ", round(i/nrow(array_grid)*100), "%")
  }



  model_time <- Sys.time()
  set.seed(123)
  m <- xgb.cv(
    data = X,
    label = Y,
    nrounds = 5000,
    objective = objective,
    metrics = "rmse",
    early_stopping_rounds = 50,
    nfold = 10,
    nthread = 1,
    verbose = 0,
    params = list(
      eta = array_grid$eta[i]
    )
  )
  array_grid$rmse[i] <- min(m$evaluation_log$test_rmse_mean)
  array_grid$trees[i] <- m$best_iteration

  total_time <- Sys.time() - model_time
  message("Model time: ")
  print(total_time)

}







n_by_response <- hyper_grid |>
  group_by(responses) |>
  count() |>
  pull(n) |>
  unique()

n_models_per_array <- 105

array_nums_1 <- rep(seq(1, ceiling(n_by_response / n_models_per_array)), each = n_models_per_array)
array_nums_2 <- array_nums_1 + max(array_nums_1)
array_nums_3 <- array_nums_1 + max(array_nums_1)*2
array_nums_4 <- array_nums_1 + max(array_nums_1)*3
max(array_nums_4)

hyper_grid <- hyper_grid |>
  mutate(array = c(array_nums_1, array_nums_2, array_nums_3, array_nums_4))

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.numeric(args[1])
if(is.na(task_id)) task_id <- 209
message("task id: ", task_id)

### start at 209!!

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

objective <- if_else(y == "mbias_density_class", "binary:logistic", "reg:squarederror")

message("Fitting xgBoost...")
start_time <- Sys.time()
cl <- makeCluster(8)
registerDoParallel(cl) # register a parallel backend
clusterExport(cl, c('X' ,'Y', 'array_grid', 'objective')) # import objects outside

out <- foreach(i = 1:3, .packages = c("xgboost"), .combine = rbind) %dopar% {

  m <- xgb.cv(
    data = X,
    label = Y,
    nrounds = 5000,
    objective = objective,
    metrics = "rmse",
    early_stopping_rounds = 50,
    nfold = 10,
    # nthread = 10,
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
  # define grid for output
  out_grid <- array_grid[i,]
  out_grid$rmse <- min(m$evaluation_log$test_rmse_mean)
  out_grid$trees <- m$best_iteration
  out_grid
}

stopCluster(cl)

total_time <- Sys.time() - start_time
message("Elapsed (foreach) time: ")
print(total_time)


# grid search
start_time <- Sys.time()
for(i in seq_len(nrow(array_grid))) {

  if(i %% 10 == 0){
    message("[", i, "/", nrow(array_grid), "] ", round(i/nrow(array_grid)*100), "%")
  }

  model_time <- Sys.time()
  set.seed(123)
  m <- xgb.cv(
    data = X,
    label = Y,
    nrounds = 5000,
    objective = objective,
    metrics = "rmse",
    early_stopping_rounds = 50,
    nfold = 10,
    nthread = 10,
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

  total_time <- Sys.time() - model_time
  message("Model time: ")
  print(total_time)

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


