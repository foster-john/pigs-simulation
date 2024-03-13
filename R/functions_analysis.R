

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
    repeats = 5
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

xgb_fit <- function(df){

  train(
    y ~ .,
    data = df,
    method = "xgbLinear",
    trControl = trainControl(method = "cv"),
    tuneGrid = hyper_grid,
    metric = "RMSE"
  )
}


my_recipe <- function(df_train, df_test){
  require(rsample)
  require(recipes)

  blueprint <- recipe(y ~ ., data = df_train) |>
    step_dummy(all_nominal_predictors()) |>
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

  } else if(ml == "xgb"){

    message("Fitting extreme gradient boosting")
    fit <- xgb_fit(df_train)

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
      pred = pred,
      data_train = df_train,
      data_test = df_test
    ),
    dest)

}

subset_rename <- function(df, y, prop = 0.7){

  require(dplyr)
  require(rsample)
  set.seed(456)

  dat <- df |>
    ungroup() |>
    mutate(simulation_id = stringr::str_extract(property_id, "[[:graph:]]*(?=-[[:digit:]]*$)"),
           methods_used = as.factor(methods_used)) |>
    rename(y = all_of(y),
           take = sum_take_density) |>
    select(-contains("density"), -contains("abundance"), -PPNum, -property, -property_id,
           -extinct, -recovered, -obs_flag, -sum_take, -contains("per_unit"))

  if(y != "mbias_density"){
    dat <- dat |> mutate(y = log(y))
  }

  simulations <- dat |> pull(simulation_id) |> unique()

  n_simulations <- length(simulations)
  n_training <- floor(n_simulations * prop)

  samps <- sample.int(n_simulations, n_training)
  training_simulations <- simulations[samps]

  train <- dat |>
    filter(simulation_id %in% training_simulations) |>
    select(-simulation_id)
  test <- dat |>
    filter(!simulation_id %in% training_simulations) |>
    select(-simulation_id)

  list(train = train, test = test)

}
