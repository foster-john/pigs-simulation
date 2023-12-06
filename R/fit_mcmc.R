fit_mcmc <- function(modelCode, data, constants, n_iter, n_chains, custom_samplers = NULL, monitors_add = NULL){

  require(nimble)

  source("R/calc_log_potential_area.R")

  Rmodel <- nimbleModel(
    code = modelCode,
    constants = constants,
    data = data,
    inits = inits(data, constants),
    calculate = TRUE
  )

  for(i in 1:constants$n_survey){
    N_model <- Rmodel$N[constants$p_property_idx[i], constants$p_pp_idx[i]]
    n <- N_model - data$y_sum[i]
    if(n < 0){
      Rmodel$N[constants$p_property_idx[i], constants$p_pp_idx[i]] <- N_model + n^2
    }
  }

  Rmodel$initializeInfo()

  # default MCMC configuration
  mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)

  if(!is.null(monitors_add)){
    mcmcConf$addMonitors(monitors_add)
  }

  if(!is.null(custom_samplers)){
    for(i in seq_len(nrow(custom_samplers))){
      node <- custom_samplers$node[i]
      type <- custom_samplers$type[i]
      mcmcConf$removeSampler(node)
      mcmcConf$addSampler(node, type)
    }
  }

  for(i in 1:5){
    node <- paste0("beta_p[", i, ", ", 1:constants$m_p, "]")
    node <- c(paste0("beta1[", i, "]"), node)
    mcmcConf$removeSampler(node)
    mcmcConf$addSampler(node, "AF_slice")
  }

  mcmcConf$printSamplers(byType = TRUE)

  Rmcmc <- buildMCMC(mcmcConf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc)

  runMCMC(
    Cmcmc,
    niter = n_iter,
    nchains = n_chains,
    nburnin = n_iter / 2,
    samplesAsCodaMCMC = TRUE
  )
}
