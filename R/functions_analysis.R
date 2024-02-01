


fit_glm_null <- function(df, outname){

  require(lme4)
  m_list <- list()

  message("  GLM - null")
  m_list$glm_0.1 <-
    glmer(
      nrmse ~ (1 | methods_used),
      family = Gamma(link = "log"),
      data = df
    )
  m_list$glm_0.2 <-
    glmer(
      nrmse ~ (1 | method),
      family = Gamma(link = "log"),
      data = df
    )
  m_list$glm_0.3 <-
    glmer(
      nrmse ~ (1 | methods_used) + (1 | method),
      family = Gamma(link = "log"),
      data = df
    )

  write_rds(m_list, outname)
  message("   Done")
  return(m_list)
}

fit_glm_method <- function(df, outname){

  require(lme4)
  m_list <- list()

  message("  GLM - null")
  m_list$glm_1.1 <-
    glmer(
      nrmse ~ (1 | method) +
        (mean_effort_method | method) +
        (mean_trap_count_method | method),
      family = Gamma(link = "log"),
      data = df
    )
  m_list$glm_1.2 <-
    glmer(
      nrmse ~ (1 | method) +
        (sum_effort_method | method) +
        (sum_trap_count_method | method),
      family = Gamma(link = "log"),
      data = df
    )

  write_rds(m_list, outname)
  message("   Done")
  return(m_list)
}

fit_gam_null <- function(df, outname){

  require(mgcv)
  m_list <- list()

  message("  GAM - null")
  gam <-
    gam(
      nrmse ~ s(methods_used, bs = "re"),
      family = Gamma(link = "log"),
      data = df)

  m_list$gam_0 <- gam
  write_rds(m_list, outname)
  message("   Done")
  return(gam)
}


fit_glm_individual <- function(df, outname){

  require(lme4)
  m_list <- list()

  message("  GLM - return_interval")
  m_list$glm_1.1 <-
    glmer(
      nrmse ~ (1 | methods_used) + return_interval,
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - med_density")
  m_list$glm_1.2 <-
    glmer(
      nrmse ~ (1 | methods_used) + med_density,
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - delta")
  m_list$glm_1.3 <-
    glmer(nrmse ~ (1 | methods_used) + delta,
          family = Gamma(link = "log"),
          data = df)
  message("   Done")

  message("  GLM - n_reps")
  m_list$glm_1.4 <-
    glmer(nrmse ~ (1 | methods_used) + n_reps,
          family = Gamma(link = "log"),
          data = df)
  message("   Done")

  message("  GLM - n_methods_used")
  m_list$glm_1.5 <-
    glmer(
      nrmse ~ (1 | methods_used) + n_methods_used,
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - sum_take_density")
  m_list$glm_1.6 <-
    glmer(
      nrmse ~ (1 | methods_used) + sum_take_density,
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - property_area")
  m_list$glm_1.7 <-
    glmer(
      nrmse ~ (1 | methods_used) + property_area,
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - mean_effort")
  m_list$glm_1.8 <-
    glmer(
      nrmse ~ (1 | methods_used) + mean_effort,
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - mean_trap_count")
  m_list$glm_1.9 <-
    glmer(
      nrmse ~ (1 | methods_used) + mean_trap_count,
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - sum_effort")
  m_list$glm_1.10 <-
    glmer(
      nrmse ~ (1 | methods_used) + sum_effort,
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - sum_trap_count")
  m_list$glm_1.11 <-
    glmer(
      nrmse ~ (1 | methods_used) + sum_trap_count,
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  write_rds(m_list, outname)
  message("GLM individual models done")
  return(m_list)

}

fit_glm_sum_take <- function(df, outname){

  require(lme4)
  m_list <- list()

  message("  GLM - sum_take_density + return_interval")
  m_list$glm_2.1 <-
    glmer(
      nrmse ~ (1 | methods_used) +
        sum_take_density +
        return_interval +
        I(sum_take_density * return_interval),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - sum_take_density + delta")
  m_list$glm_2.2 <-
    glmer(nrmse ~ (1 | methods_used) +
            sum_take_density +
            delta +
            I(sum_take_density * delta),
          family = Gamma(link = "log"),
          data = df)
  message("   Done")

  message("  GLM - sum_take_density + n_reps")
  m_list$glm_2.3 <-
    glmer(nrmse ~ (1 | methods_used) +
            sum_take_density +
            n_reps +
            I(sum_take_density * n_reps),
          family = Gamma(link = "log"),
          data = df)
  message("   Done")

  message("  GLM - sum_take_density + n_methods_used")
  m_list$glm_2.4 <-
    glmer(
      nrmse ~ (1 | methods_used) +
        sum_take_density +
        n_methods_used +
        I(sum_take_density * n_methods_used),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - sum_take_density + property_area")
  m_list$glm_2.5 <-
    glmer(
      nrmse ~ (1 | methods_used) +
        sum_take_density +
        property_area +
        I(sum_take_density * property_area),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - sum_take_density + sum_effort")
  m_list$glm_2.6 <-
    glmer(
      nrmse ~ (1 | methods_used) +
        sum_take_density +
        sum_effort +
        I(sum_take_density * sum_effort),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - sum_take_density + sum_trap_count")
  m_list$glm_2.7 <-
    glmer(
      nrmse ~ (1 | methods_used) +
        sum_take_density +
        sum_trap_count +
        I(sum_take_density * sum_trap_count),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  write_rds(m_list, outname)
  message("GLM sum take + interactions models done")
  return(m_list)

}

fit_glm_sum_take_area <- function(df, outname){

  require(lme4)
  m_list <- list()

  message("  GLM - sum_take_density + property area + return_interval")
  m_list$glm_3.1 <-
    glmer(
      nrmse ~ (1 | methods_used) +
        sum_take_density +
        property_area +
        return_interval +
        I(sum_take_density * property_area) +
        I(sum_take_density * return_interval) +
        I(property_area * return_interval),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - sum_take_density + property area + delta")
  m_list$glm_3.2 <-
    glmer(
      nrmse ~ (1 | methods_used) +
        sum_take_density +
        property_area +
        delta +
        I(sum_take_density * property_area) +
        I(sum_take_density * delta) +
        I(property_area * delta),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - sum_take_density + property area + n_reps")
  m_list$glm_3.3 <-
    glmer(
      nrmse ~ (1 | methods_used) +
        sum_take_density +
        property_area +
        n_reps +
        I(sum_take_density * property_area) +
        I(sum_take_density * n_reps) +
        I(property_area * n_reps),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - sum_take_density + property area + sum_effort")
  m_list$glm_3.4 <-
    glmer(
      nrmse ~ (1 | methods_used) +
        sum_take_density +
        property_area +
        sum_effort +
        I(sum_take_density * property_area) +
        I(sum_take_density * sum_effort) +
        I(property_area * sum_effort),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  message("  GLM - sum_take_density + property area + sum_trap_count")
  m_list$glm_3.5 <-
    glmer(
      nrmse ~ (1 | methods_used) +
        sum_take_density +
        property_area +
        sum_trap_count +
        I(sum_take_density * property_area) +
        I(sum_take_density * sum_trap_count) +
        I(property_area * sum_trap_count),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  write_rds(m_list, outname)
  message("GLM sum take + interactions models done")
  return(m_list)

}

fit_glm_all <- function(df, outname){

  require(lme4)
  m_list <- list()

  message("  GLM - ALL")
  m_list$glm_5 <-
    glmer(
      nrmse ~ (1 | methods_used) +
        sum_take_density +
        med_density +
        property_area +
        return_interval +
        delta +
        n_reps +
        sum_trap_count +
        sum_effort +
        n_methods_used +
        I(sum_take_density * med_density) +
        I(sum_take_density * property_area) +
        I(sum_take_density * return_interval) +
        I(sum_take_density * delta) +
        I(sum_take_density * n_reps) +
        I(sum_take_density * sum_trap_count) +
        I(sum_take_density * sum_effort) +
        I(sum_take_density * n_methods_used) +
        I(med_density * property_area) +
        I(med_density * return_interval) +
        I(med_density * delta) +
        I(med_density * n_reps) +
        I(med_density * sum_trap_count) +
        I(med_density * sum_effort) +
        I(med_density * n_methods_used) +
        I(property_area * return_interval) +
        I(property_area * delta) +
        I(property_area * n_reps) +
        I(property_area * sum_trap_count) +
        I(property_area * sum_effort) +
        I(property_area * n_methods_used) +
        I(return_interval * delta) +
        I(return_interval * n_reps) +
        I(return_interval * sum_trap_count) +
        I(return_interval * sum_effort) +
        I(return_interval * n_methods_used) +
        I(delta * n_reps) +
        I(delta * sum_trap_count) +
        I(delta * sum_effort) +
        I(delta * n_methods_used) +
        I(n_reps * sum_trap_count) +
        I(n_reps * sum_effort) +
        I(n_reps * n_methods_used) +
        I(sum_trap_count * sum_effort) +
        I(sum_trap_count * n_methods_used) +
        I(sum_effort * n_methods_used),
      family = Gamma(link = "log"),
      data = df
    )
  write_rds(m_list, outname)
  message("  Done")
  return(m_list)
}

fit_gam_individual <- function(df, outname){

  require(mgcv)
  m_list <- list()

  message("  GAM - return_interval")
  m_list$gam_1.1 <-
    gam(
      nrmse ~ s(methods_used, bs = "re") + s(return_interval, bs = "cr"),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")
  message("  GAM - med_density")
  m_list$gam_1.2 <-
    gam(
      nrmse ~ s(methods_used, bs = "re") + s(med_density, bs = "cr"),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")
  message("  GAM - delta")
  m_list$gam_1.3 <-
    gam(
      nrmse ~ s(methods_used, bs = "re") + s(delta, bs = "cr"),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")
  message("  GAM - n_reps")
  m_list$gam_1.4 <-
    gam(
      nrmse ~ s(methods_used, bs = "re") + s(n_reps, bs = "cr"),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")
  message("  GAM - sum_take_density")
  m_list$gam_1.5 <-
    gam(
      nrmse ~ s(methods_used, bs = "re") + s(sum_take_density, bs = "cr"),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")
  message("  GAM - property_area")
  m_list$gam_1.6 <-
    gam(
      nrmse ~ s(methods_used, bs = "re") + s(property_area, bs = "cr"),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")
  message("  GAM - mean_effort")
  m_list$gam_1.7 <-
    gam(
      nrmse ~ s(methods_used, bs = "re") + s(mean_effort, bs = "cr"),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")
  message("  GAM - mean_trap_count")
  m_list$gam_1.8 <-
    gam(
      nrmse ~ s(methods_used, bs = "re") + s(mean_trap_count, bs = "cr"),
      family = Gamma(link = "log"),
      data = df
    )
  message("   Done")

  write_rds(m_list, outname)
  message("GAM individual models done")
  return(m_list)
}


# gam <- gam(nrmse ~
#              s(methods_used, bs = "re") +
#              s(return_interval, bs = "cr") +
#              s(med_density, bs = "cr") +
#              s(delta, bs = "cr") +
#              s(n_reps, bs = "cr") +
#              s(n_methods_used, bs = "cr") +
#              s(sum_take_density, bs = "cr") +
#              s(property_area, bs = "cr") +
#              s(mean_effort, bs = "cr") +
#              s(mean_trap_count, bs = "cr"),
#            family = Gamma(link = "log"),
#            data = df))
