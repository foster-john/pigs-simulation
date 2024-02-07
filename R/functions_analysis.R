
f_nrmse <- formula(
  y ~ (1 | methods_used) +
    property_area +
    total_take_density +
    delta +
    # effort +                              # removed after 1st round
    # unit_count +                          # removed after 1st round
    n_reps_pp +
    I_property_area_x_total_take_density +
    # I_property_area_x_delta +             # removed after 1st round
    I_property_area_x_unit_count +
    I_property_area_x_n_reps_pp +
    # I_property_area_x_effort +            # removed after 2nd round
    I_total_take_density_x_delta +
    # I_total_take_density_x_unit_count +   # removed after 1st round
    I_total_take_density_x_n_reps_pp +
    I_total_take_density_x_effort +
    # I_delta_x_unit_count +                # removed after 1st round
    I_delta_x_n_reps_pp +
    I_delta_x_effort
    # I_unit_count_x_n_reps_pp +            # removed after 1st round
    # I_unit_count_x_effort +               # removed after 1st round
    # I_n_reps_pp_x_effort                  # removed after 1st round
)

f_mpe <- formula(
  y ~ (1 | methods_used) +
    property_area +
    total_take_density +
    # delta +                               # removed after 1st round
    effort +
    # unit_count +                          # removed after 1st round
    n_reps_pp +
    I_property_area_x_total_take_density +
    # I_property_area_x_delta +             # removed after 1st round
    I_property_area_x_unit_count +
    I_property_area_x_n_reps_pp +
    # I_property_area_x_effort +            # removed after 1st round
    I_total_take_density_x_delta +
    # I_total_take_density_x_unit_count +   # removed after 1st round
    I_total_take_density_x_n_reps_pp +
    I_total_take_density_x_effort +
    # I_delta_x_unit_count +                # removed after 1st round
    I_delta_x_n_reps_pp +
    I_delta_x_effort
    # I_unit_count_x_n_reps_pp +            # removed after 1st round
    # I_unit_count_x_effort +               # removed after 1st round
    # I_n_reps_pp_x_effort                  # removed after 1st round
)

f_bias <- formula(
  y ~ (1 | methods_used) +
    property_area +
    total_take_density +
    delta +
    # effort +                              # removed after 1st round
    # unit_count +                          # removed after 1st round
    # n_reps_pp +                           # removed after 1st round
    I_property_area_x_total_take_density +
    I_property_area_x_delta +
    # I_property_area_x_unit_count +        # removed after 1st round
    I_property_area_x_n_reps_pp +
    # I_property_area_x_effort +            # removed after 1st round
    I_total_take_density_x_delta +
    I_total_take_density_x_unit_count +
    I_total_take_density_x_n_reps_pp +
    # I_total_take_density_x_effort +       # removed after 1st round
    # I_delta_x_unit_count +                # removed after 1st round
    # I_delta_x_n_reps_pp +                 # removed after 1st round
    I_delta_x_effort
    # I_unit_count_x_n_reps_pp +            # removed after 1st round
    # I_unit_count_x_effort +               # removed after 1st round
    # I_n_reps_pp_x_effort                  # removed after 1st round
)

f <- formula(
  y ~ methods_used +
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

f_gam <- formula(
  y ~ s(methods_used, bs = "re") +
    s(property_area) +
    s(total_take_density) +
    s(delta) +
    s(effort) +
    s(unit_count) +
    s(n_reps_pp) +
    s(I_property_area_x_total_take_density) +
    s(I_property_area_x_delta) +
    s(I_property_area_x_unit_count) +
    s(I_property_area_x_n_reps_pp) +
    s(I_property_area_x_effort) +
    s(I_total_take_density_x_delta) +
    s(I_total_take_density_x_unit_count) +
    s(I_total_take_density_x_n_reps_pp) +
    s(I_total_take_density_x_effort) +
    s(I_delta_x_unit_count) +
    s(I_delta_x_n_reps_pp) +
    s(I_delta_x_effort) +
    s(I_unit_count_x_n_reps_pp) +
    s(I_unit_count_x_effort) +
    s(I_n_reps_pp_x_effort)
)

fit_glm_all <- function(df, y, vars, path){

  require(lme4)
  require(dplyr)
  require(mgcv)

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

  if(y != "bias"){
    data <- data |> mutate(y = log(y))
  }

  n <- gam(y ~ s(methods_used, bs = "re"), data = data)
  fit <- gam(f_gam, data = data)

  filename <- paste(y, vars, sep = "-")

  outname <- file.path(path, paste0(filename, "-null.rds"))
  write_rds(list(fit = n, data = data), outname)

  outname <- file.path(path, paste0(filename, ".rds"))
  write_rds(list(fit = fit, data = data), outname)

  message("  Inital fit done")
  message("  Fit warnings:")
  print(warnings())

  # message(" Step backward")
  # step_back <- step(fit, direction = "backward")
  #
  # outname <- file.path(path, paste0(filename, "-stepback.rds"))
  # write_rds(step_back, outname)

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

