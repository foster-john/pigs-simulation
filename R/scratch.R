library(nimble)
library(parallel)
library(coda)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

path <- "out/hpc/density_1.4"
task <- "20"


rds <- read_rds(file.path(path, task, "simulation_data.rds"))
rds$psrf

burn <- "parameters_burnin.rds"
params <- read_rds(file.path(path, task, burn))
params <- params[[1]]

plot(params[,"psi_phi"])

rds$posterior_samples |>
  pull("psi_phi") |>
  # exp() |>
  hist(breaks = 50)

rds$N$property |> max()

phi_mu <- tibble(
  phi = rds$posterior_samples |> pull("phi_mu"),
  psi = rds$posterior_samples |> pull("psi_phi")
) |>
  mutate(a = phi * psi,
         b = (1-phi) * psi,
         phi_mu = rbeta(n(), a, b))

hist(phi_mu$phi_mu, breaks = 50)
quantile(phi_mu$phi_mu, c(0.025, 0.5, 0.975))


rds$N |>
  mutate(property = as.character(property)) |>
  ggplot() +
  aes(x = PPNum, y = N, color = property) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

str(rds)


n_samps <- rds$posterior_samples |>
  select(contains("xn"))

n_median <- apply(n_samps, 2, median)
n_low <- apply(n_samps, 2, quantile, 0.025)
n_high <- apply(n_samps, 2, quantile, 0.975)

tibble(
  observed = rds$N$N,
  predicted = n_median,
  ymin = n_low,
  ymax = n_high
) |>
  ggplot() +
  aes(x = observed, y = predicted, ymin = ymin, ymax = ymax) +
  geom_point() +
  # geom_linerange() +
  geom_abline(intercept = 0, slope = 1) +
  # coord_cartesian(ylim = c(0, 250)) +
  theme_bw()

