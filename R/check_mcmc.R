check_mcmc <- function(samples, nodes_check, n_mcmc, dest){

  require(coda)

  all_nodes <- colnames(samples[[1]])
  n_iter <- nrow(samples[[1]])

  j <- unlist(lapply(nodes_check, function(x) grep(x, all_nodes)))

  params <- samples[,j]

  message("Calculating PSRF...")
  psrf <- gelman.diag(params, multivariate = FALSE)
  print(psrf)

  message("Calculating effective samples...")
  effective_samples <- effectiveSize(params)
  print(effective_samples)

 # message("Calculating burnin...")
 # ff <- tempfile()
 # png(filename = ff)
 # GBR <- gelman.plot(params)
 # dev.off()
 # unlink(ff)

 # burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
 # print(burnin)

 # if(is.na(burnin) | burnin >= 0.9*nrow(params[[1]])){
    burnin <- round(n_iter / 2)
 # }

  samples_burn_mcmc <- window(samples, start = burnin)
  params_burn <- samples_burn_mcmc[, j]

 # message("Creating traceplots...")
 # png(filename = file.path(dest, "mcmcTimeseries%03d.png"))
 # plot(params_burn)
 # dev.off()
 # message("  done")

  samples_burn_mat <- as.matrix(samples_burn_mcmc)
  draws <- sample.int(nrow(samples_burn_mat), n_mcmc, replace = TRUE)

  samples_draw <- as.data.frame(samples_burn_mat[draws, ])

  list(
    posterior_samples = samples_draw,
    psrf = psrf$psrf,
    effective_samples = effective_samples,
    burnin = burnin,
    converged = max(psrf$psrf[, "Upper C.I."]) <= 1.1,
    bad_mcmc = any(is.na(psrf$psrf)) | any(psrf$psrf > 3)
  )

}
