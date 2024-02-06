
f_nmrmse <- formula(
  y ~ (1 | methods_used) +
    property_area +
    total_take_density +
    delta +
    unit_count +
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
    property_area +
    total_take_density +
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
    I(total_take_density * delta) +
    I(total_take_density * unit_count) +
    I(total_take_density * n_reps_pp) +
    I(delta * effort)
)

fit_glm_all <- function(df, y, effort, agg, path){

  require(lme4)
  require(dplyr)

  subset_rename <- function(df, y, effort, agg){
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

    agg_select <- if_else(agg == "mean", "sum", "mean")
    effort_select <- if_else(effort == "per_unit", "raw", "per_unit")

    df |>
      select(-contains(agg_select), -contains(effort_select)) |>
      rename(effort = contains("effort"),
             unit_count = contains("unit_count"))
  }

  data <- subset_rename(df, y, effort, agg)

  message("y = ", y)
  message("effort = ", effort)
  message("agg = ", agg)

  # data |>
  #   select(effort, unit_count, n_reps_pp) |>
  #   cor() |>
  #   print()

  if(y == "bias") {
    n <- lmer(y ~ (1 | methods_used), data = data)
    fit <- lmer(f_bias, data = data)
  } else if(y == "nrmse") {
    n <- glmer(y ~ (1 | methods_used) , family = Gamma(link = "log"), data = data)
    fit <- glmer(f_nmrmse , family = Gamma(link = "log"), data = data)
  } else if(y == "mpe") {
    n <- glmer(y ~ (1 | methods_used) , family = Gamma(link = "log"), data = data)
    fit <- glmer(f_mpe , family = Gamma(link = "log"), data = data)
  }

  filename <- paste(y, effort, agg, sep = "-")

  outname <- file.path(path, paste0(filename, "-null.rds"))
  write_rds(list(fit = n, data = data), outname)

  outname <- file.path(path, paste0(filename, "-cut1.rds"))
  write_rds(list(fit = fit, data = data), outname)

  message("  Inital fit done")
  message("  Fit warnings:")
  print(warnings())

  # message("Dredge")
  # oop <- options(na.action = "na.fail")
  # dd <- MuMIn::dredge(fit)
  #
  # outname <- file.path(path, paste0(filename, "-dredge.rds"))
  # write_rds(dd, outname)
  #
  # message("  Dredge done")
  # message("  Dredge warnings:")
  # print(warnings())

  return(list(fit = fit, dredge = dd))
}

