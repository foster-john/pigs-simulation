

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
    select(-contains("density"), -contains("abundance"), -PPNum, -property,
           -extinct, -recovered, -obs_flag)

  if(y %in% c("nm_rmse_density", "mpe_density", "med_density", "var_density")){
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

  # create strata by decile
  # each property will belong to a decile of each land cover variable
  df_strata <- dat |>
    mutate(canopy_strata = make_strata(c_canopy, breaks = 10),
           rugged_strata = make_strata(c_rugged, breaks = 10),
           road_den_strata = make_strata(c_road_den, breaks = 10)) |>
    select(property_id, contains("strata")) |>
    distinct()

  col_sample <- function(dfs, col){
    min_sample <- dfs |>
      pull(all_of(col)) |>
      table() |>
      min()

    dfs |>
      group_by(.data[[col]]) |>
      slice_sample(n = min(min_sample, 15)) |>
      ungroup() |>
      pull(property_id)
  }

  canopy_sample <- col_sample(df_strata, "canopy_strata")
  rugged_sample <- col_sample(df_strata, "rugged_strata")
  road_den_sample <- col_sample(df_strata, "road_den_strata")

  props <- c(canopy_sample, rugged_sample, road_den_sample)
  df_sample <- df_strata |>
    filter(property_id %in% unique(props))

  message("Number of properties in each canopy strata:")
  print(table(df_sample$canopy_strata))
  message("Number of properties in each rugged strata:")
  print(table(df_sample$rugged_strata))
  message("Number of properties in each road density strata:")
  print(table(df_sample$road_den_strata))

  train <- dat |>
    filter(property_id %in% props)

  training_simulations <- unique(train$simulation_id)

  test <- dat |>
    filter(!property_id %in% props,
           !simulation_id %in% training_simulations)

  testing_simulations <- unique(test$simulation_id)

  message("Number of simulations in training data: ", length(training_simulations))
  message("Number of simulations in testing data: ", length(testing_simulations))

  total_simulations <- length(unique(dat$simulation_id))
  testthat::expect_equal(total_simulations,
                         length(training_simulations) + length(testing_simulations))

  list(train = train |> select(-simulation_id, -property_id),
       test = test |> select(-simulation_id, -property_id))

}
