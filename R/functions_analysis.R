

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

subset_rename <- function(df, y, prop = 0.8){

  require(dplyr)
  require(rsample)
  set.seed(456)

  dat <- df |>
    ungroup() |>
    mutate(simulation_id = stringr::str_extract(property_id, "[[:graph:]]*(?=-[[:digit:]]*$)"),
           methods_used = as.factor(methods_used)) |>
    rename(y = all_of(y),
           sum_take_d = sum_take_density) |>
    group_by(property_id) |>
    mutate(delta = c(NA, diff(PPNum))) |>
    ungroup() |>
    select(-contains("density"), -contains("abundance"), -PPNum, -property, -property_id,
           -extinct, -recovered, -obs_flag)

  if(y %in% c("nm_rmse_density", "mpe_density")){
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
