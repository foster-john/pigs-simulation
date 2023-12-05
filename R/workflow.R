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



