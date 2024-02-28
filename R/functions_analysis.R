

ranger_fit <- function(df){

  require(ranger)
  require(caret)

  # number of features
  n_features <- length(setdiff(names(df), "y"))

  # create hyperparameter grid
  hyper_grid <- expand.grid(
    splitrule = "variance",
    mtry = floor(n_features * c(.15, .25, .333, .4)),
    min.node.size = c(1, 3, 5, 10, 15)
    # replace = c(TRUE, FALSE),
    # sample.fraction = c(.5, .63, .8),
    # rmse = NA
  )

  train(
    y ~ .,
    data = df,
    method = "ranger",
    trControl = trainControl(method = "cv"),
    tuneGrid = hyper_grid,
    metric = "RMSE",
    importance = "impurity"
  )

}

knn_fit <- function(df){

  require(caret)

  cv <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )


  # Create a hyperparameter grid search
  hyper_grid <- expand.grid(
    k = floor(seq(1, nrow(df)/3, length.out = 20))
  )

  # Fit knn model and perform grid search
  train(
    y ~ .,
    data = df,
    method = "knn",
    trControl = cv,
    tuneGrid = hyper_grid,
    metric = "RMSE"
  )
}

my_recipe <- function(df){
  require(rsample)
  require(recipes)

  split <- initial_split(df, prop = 0.6,
                         strata = "methods_used")

  df_train <- training(split)
  df_test <- testing(split)

  blueprint <- recipe(y ~ ., data = df_train) |>
    step_dummy(all_nominal_predictors()) |>
    step_interact(terms = ~ starts_with("methods_used"):all_numeric_predictors()) |>
    step_interact(terms = ~ property_area:take) |>
    step_interact(terms = ~ property_area:delta) |>
    step_interact(terms = ~ property_area:mean_effort) |>
    step_interact(terms = ~ property_area:sum_effort) |>
    step_interact(terms = ~ property_area:mean_unit_count) |>
    step_interact(terms = ~ property_area:sum_unit_count) |>
    step_interact(terms = ~ property_area:n_reps_pp) |>
    step_interact(terms = ~ take:delta) |>
    step_interact(terms = ~ take:mean_effort) |>
    step_interact(terms = ~ take:sum_effort) |>
    step_interact(terms = ~ take:mean_unit_count) |>
    step_interact(terms = ~ take:sum_unit_count) |>
    step_interact(terms = ~ take:n_reps_pp) |>
    step_interact(terms = ~ delta:mean_effort) |>
    step_interact(terms = ~ delta:sum_effort) |>
    step_interact(terms = ~ delta:mean_unit_count) |>
    step_interact(terms = ~ delta:sum_unit_count) |>
    step_interact(terms = ~ delta:n_reps_pp) |>
    step_interact(terms = ~ n_reps_pp:mean_effort) |>
    step_interact(terms = ~ n_reps_pp:sum_effort) |>
    step_interact(terms = ~ n_reps_pp:mean_unit_count) |>
    step_interact(terms = ~ n_reps_pp:sum_unit_count) |>
    step_interact(terms = ~ mean_effort:mean_unit_count) |>
    step_interact(terms = ~ mean_effort:sum_unit_count) |>
    step_interact(terms = ~ sum_effort:mean_unit_count) |>
    step_interact(terms = ~ sum_effort:sum_unit_count) |>
    step_nzv(all_predictors()) |>
    step_center(all_numeric_predictors()) |>
    step_scale(all_numeric_predictors())

  prepare <- prep(blueprint, training = df_train)

  baked_train <- bake(prepare, new_data = df_train)
  baked_test <- bake(prepare, new_data = df_test)

  return(list(df_train = baked_train, df_test = baked_test))

}

fit_ml <- function(df, ml, dest){

  require(dplyr)
  require(doParallel)

  df_ls <- my_recipe(df)
  df_train <- df_ls$df_train
  df_test <- df_ls$df_test

  message("Make cluster")
  cl <- makePSOCKcluster(5)
  registerDoParallel(cl)

  if(ml == "ranger"){

    message("Fitting random forest...")
    fit <- ranger_fit(df_train)

  } else if(ml == "knn"){

    message("Fitting k-nearest neighbors")
    fit <- knn_fit(df_train)

  }

  pred <- predict(fit, df_test)
  stopCluster(cl)

  message(paste0("   ", ml, " fit done"))
  message("  Fit warnings:")
  print(warnings())

  message("==== Fit ====")
  print(fit)

  message("==== Prediction ====")
  print(postResample(pred, df_test$y))


  write_rds(
    list(
      fit = fit,
      data_train = df_train,
      data_test = df_test
    ),
    dest)

}

subset_rename <- function(df, y, outlier_quant = 0.995){

  require(dplyr)

  outlier <- df |>
    pull(all_of(y)) |>
    quantile(outlier_quant)

  dat <- df |>
    ungroup() |>
    rename(y = all_of(y),
           take = sum_take_density) |>
    # filter(density > 0,
    #        y < outlier) |>
    select(-contains("density"), -contains("abundance"), -PPNum, -property, -property_id,
           -extinct, -recovered, -obs_flag, -sum_take, -contains("per_unit"))

  if(y != "mbias_density"){
    dat <- dat |> mutate(y = log(y))
  }

  return(dat)

}
