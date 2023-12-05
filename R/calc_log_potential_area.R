calc_log_potential_area <- nimbleFunction(
  run = function(
    log_rho = double(1),
    log_gamma = double(1),
    p_unique = double(1),
    log_effort_per = double(0),
    effort_per = double(0),
    n_trap_m1 = double(0),
    log_pi = double(0),
    method = double(0)
  ){
    m <- method

    if(m == 1){ # firearms
      log_potential_area <- log_rho[m] +
        log_effort_per -
        log(1 + (p_unique[m] * n_trap_m1))
    } else if(m == 2 | m == 3){ # fixed wing and helicopter
      log_potential_area <- log_rho[m] + log_effort_per
    } else if(m == 4 | m == 5){
      log_potential_area <- log_pi +
        (2 * (log_rho[m] + log_effort_per -
                log(exp(log_gamma[m-3]) + effort_per))) +
        log(1 + (p_unique[m-2] * n_trap_m1))
    }
    return(log_potential_area)
    returnType(double(0))
  },
  buildDerivs = TRUE
)
