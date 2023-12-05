# --------------------------------------------------------------------
#
# Workflow script for pig removal simulations
#
# John Foster
#
# General workflow:
#
# 1. Simulate/bootstrap data
#   - Load MIS data
#   - 1-method properties
#   - 2- to 5-method properties
# 2. Simulate eco/take dynamics
# 3. Fit MCMC
# 4. Check MCMC
# 5. Write summary statistics
#
# --------------------------------------------------------------------

library(tidyverse)
library(nimble)
library(config)

config <- get(config = "default")

# -----------------------------------------------------------------
# Load MIS data ----
# -----------------------------------------------------------------

df <- read_rds("../pigs-statistical/data/insitu/MIS_4weekPP.rds")
df <- df |>
  select(-method) |>
  rename(property = agrp_prp_id,
         method = Method)

# need a lookup table for property IDs and how may methods they employ ----
n_method_lookup <- df |>
  select(property, method) |>
  distinct() |>
  count(property)

# get the proportion of properties that use n methods
# for determining the number of properties to simulate
n_rel <- n_method_lookup |>
  count(n, name = "n_sum") |>
  mutate(rel_prop = n_sum / sum(n_sum),
         n_simulate = ceiling(rel_prop * 100))


n_pp <- config$n_pp                # the number of primary periods to simulate
n_prop <- sum(n_rel$n_simulate)    # total number of properties to simulate

# -----------------------------------------------------------------
# 1-method properties ----
# -----------------------------------------------------------------

source("R/one_method_properties.R")
method_1 <- one_method_properties(df, n_rel$n_simulate[1], n_pp)

# -----------------------------------------------------------------
# n-method properties ----
# -----------------------------------------------------------------

source("R/n_method_properties.R")

## use map here?
method_2 <- n_method_properties(df, n_rel$n_simulate[2], 2, n_pp)
method_3 <- n_method_properties(df, n_rel$n_simulate[3], 3, n_pp)
method_4 <- n_method_properties(df, n_rel$n_simulate[4], 4, n_pp)
method_5 <- n_method_properties(df, n_rel$n_simulate[5], 5, n_pp)

properties <- c(
  method_1,
  method_2,
  method_3,
  method_4,
  method_5
)

n_properties <- length(properties)
assertthat::are_equal(n_properties, sum(n_rel$n_simulate))

# -----------------------------------------------------------------
# Simulate eco/take dynamics ----
# -----------------------------------------------------------------

## the first thing we need to do is assign properties to counties
max_counties <- 15
weights <- runif(max_counties)
norm_weights <- weights / sum(weights)
properties_per_county <- rmulti(1, n_properties, norm_weights)
properties_per_county <- properties_per_county[properties_per_county > 0]
n_county <- length(properties_per_county)

## randomly assign which county each property goes in
order_county <- tibble(
  order_property = sample.int(n_properties),
  order_county = rep(1:n_county, times = properties_per_county)
) |>
  arrange(order_property) |>
  pull(order_county)

## we need covariate data
land_cover <- matrix(rnorm(n_county * 3), n_county, 3)

## we need the effect of X on detection probability p
beta_p <- matrix(rnorm(20), 5, 4)

# method lookup table (data model parameters)
method_lookup <- tibble(
  idx = 1:5,
  method = c("Firearms", "Fixed wing", "Helicopter", "Snares", "Traps"),
  p_unique = c(runif(1), 0, 0, runif(2)),
  rho = c(runif(1, 0.01, 5), # firearms; p_mu[1]
          runif(1, 1, 30),   # fixed wing
          runif(1, 1, 30),   # helicopter
          runif(1, 0.01, 5),  # snare; gamma[1], p_mu[2]
          runif(1, 0.01, 5)), # traps; gamma[2], p_mu[3]
  gamma = c(0, 0, 0, rgamma(1, 7.704547, 4.41925), rgamma(1, 3.613148, 3.507449))
)


phi_mu <- config$phi_mu
psi_phi <- config$psi_phi
start_density <- config$start_density
source("R/eco_dynamics.R")

# x <- 1
# simulate_dm(
#   properties[[x]],
#   order_county[x],
#   phi_mu,
#   psi_phi,
#   land_cover[order_county[x], ],
#   beta_p,
#   start_density,
#   method_lookup
# )

all_dynamics <- 1:n_properties |>
  map(
    \(x) simulate_dm(
      properties[[x]],
      x,
      order_county[x],
      phi_mu,
      psi_phi,
      land_cover[order_county[x], ],
      beta_p,
      start_density,
      method_lookup
    )
  ) |>
  list_rbind()

N <- all_dynamics |>
  select(PPNum, N, property, county, property_area) |>
  distinct() |>
  mutate(density = N / property_area)

take <- all_dynamics |>
  filter(!is.na(take))


# -----------------------------------------------------------------
# Fit MCMC ----
# -----------------------------------------------------------------

source("R/prep_nimble.R")
nimble_data <- prep_nimble(N, take, land_cover)
constants <- nimble_data$constants
data <- nimble_data$data

source("R/inits.R")
inits_test <- inits(data, constants)

custom_samplers <- tribble(
  ~node,            ~type,
  "log_mean_ls",    "slice",
  "phi_mu",         "slice",
  "psi_phi",        "slice",
  "log_rho",        "AF_slice"
)

source("R/calc_log_potential_area.R")
source("R/model_removal_dm.R")
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

monitors_add <- c("xn", "p", "log_theta")
params_check <- c(
  "beta_p",
  "beta1",
  "log_gamma",
  "log_rho",
  "phi_mu",
  "psi_phi",
  "log_mean_ls",
  "p_mu"
)

mcmcConf$addMonitors(monitors_add)

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

n_iter <- 5000
n_chains <- 3

samples <- runMCMC(
  Cmcmc,
  niter = n_iter,
  nchains = n_chains,
  nburnin = n_iter / 2,
  samplesAsCodaMCMC = TRUE
)


