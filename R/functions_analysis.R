
f_nmrmse <- formula(
  y ~ (1 | methods_used) +
    property_area +
    total_take_density +
    delta +
    n_reps_pp +
    effort +
    I(property_area * total_take_density) +
    I(property_area * delta) +
    I(property_area * unit_count) +
    I(property_area * n_reps_pp) +
    I(total_take_density * delta) +
    I(total_take_density * unit_count) +
    I(total_take_density * n_reps_pp) +
    I(total_take_density * effort) +
    I(delta * n_reps_pp) +
    I(delta * effort)
)

f_mpe <- formula(
  y ~ (1 | methods_used) +
    unit_count +
    n_reps_pp +
    effort +
    I(property_area * total_take_density) +
    I(property_area * unit_count) +
    I(property_area * n_reps_pp) +
    I(total_take_density * delta) +
    I(total_take_density * unit_count) +
    I(total_take_density * n_reps_pp) +
    I(total_take_density * effort) +
    I(delta * unit_count) +
    I(delta * n_reps_pp) +
    I(delta * effort)
)

f_bias <- formula(
  y ~ (1 | methods_used) +
    property_area +
    total_take_density +
    delta +
    n_reps_pp +
    I(property_area * total_take_density) +
    I(property_area * delta) +
    I(property_area * n_reps_pp) +
    I(total_take_density * unit_count) +
    I(total_take_density * n_reps_pp) +
    I(delta * effort)
)

f <- formula(
  y ~ (1 | methods_used) +
    property_area +
    total_take_density +
    delta +
    effort +
    unit_count +
    n_reps_pp +
    I_property_area_x_total_take_density +
    I_property_area_x_delta +
    I_property_area_x_unit_count +
    I_property_area_x_n_reps_pp +
    I_property_area_x_effort +
    I_total_take_density_x_delta +
    I_total_take_density_x_unit_count +
    I_total_take_density_x_n_reps_pp +
    I_total_take_density_x_effort +
    I_delta_x_unit_count +
    I_delta_x_n_reps_pp +
    I_delta_x_effort +
    I_unit_count_x_n_reps_pp +
    I_unit_count_x_effort +
    I_n_reps_pp_x_effort
)

fit_glm_all <- function(df, y, vars, path){

  require(lme4)
  require(dplyr)

  subset_rename <- function(df, y){
    if(y == "nrmse"){
      df <- df |>
        rename(y = nm_rmse_density) |>
        select(-mbias_density, -mpe_density)
    } else if (y == "bias"){
      df <- df |>
        rename(y = mbias_density) |>
        select(-nm_rmse_density, -mpe_density)
    } else if(y == "mpe"){
      df <- df |>
        rename(y = mpe_density) |>
        select(-mbias_density, -nm_rmse_density)
    }

  }

  data <- subset_rename(df, y)

  message("y = ", y)


  if(y == "bias") {

    n <- lmer(y ~ (1 | methods_used), data = data)
    fit <- lmer(f, data = data)

  } else {

    data <- data |> mutate(y = log(y))

    if(y == "nrmse") {

      n <- lmer(y ~ (1 | methods_used) , data = data)
      fit <- lmer(f , data = data)

    } else if(y == "mpe") {

      n <- lmer(y ~ (1 | methods_used) , data = data)
      fit <- lmer(f , data = data)

    }

  }

  filename <- paste(y, vars, sep = "-")

  outname <- file.path(path, paste0(filename, "-null.rds"))
  write_rds(list(fit = n, data = data), outname)

  outname <- file.path(path, paste0(filename, ".rds"))
  write_rds(list(fit = fit, data = data), outname)

  message("  Inital fit done")
  message("  Fit warnings:")
  print(warnings())

  message("Dredge")
  oop <- options(na.action = "na.fail")
  dd <- MuMIn::dredge(fit)

  outname <- file.path(path, paste0(filename, "-dredge.rds"))
  write_rds(dd, outname)

  message("  Dredge done")
  message("  Dredge warnings:")
  print(warnings())

  return(list(fit = fit, dredge = dd))
}

