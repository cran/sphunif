## ----setup, include = FALSE---------------------------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
  collapse = TRUE,
  message = FALSE,
  comment = "#>"
)

## ----cir_data-----------------------------------------------------------------
set.seed(987202226)
cir_data <- runif(n = 30, min = 0, max = 2 * pi)

## ----unif_test_cir------------------------------------------------------------
library(sphunif)
unif_test(data = cir_data, type = "Watson") # An htest object

## ----asymp--------------------------------------------------------------------
unif_test(data = cir_data, type = "Watson", p_value = "MC") # Monte Carlo
unif_test(data = cir_data, type = "Watson", p_value = "asymp") # Asymp. distr.

## ----avail_cir----------------------------------------------------------------
avail_cir_tests

## ----unif_test_avail_cir------------------------------------------------------
# A *list* of htest objects
unif_test(data = cir_data, type = c("PAD", "Watson"))

## ----sph_data-----------------------------------------------------------------
# Sample data on S^2
set.seed(987202226)
theta <- runif(n = 30, min = 0, max = 2 * pi)
phi <- runif(n = 30, min = 0, max = pi)
sph_data <- cbind(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi))

## ----avail_sph----------------------------------------------------------------
avail_sph_tests

## ----unif_test_avail_sph------------------------------------------------------
unif_test(data = sph_data, type = "all", p_value = "MC", M = 1e2)
unif_test(data = sph_data, type = "Rayleigh", p_value = "asymp")

## ----hyp_data-----------------------------------------------------------------
# Sample data on S^9
hyp_data <- r_unif_sph(n = 30, p = 10)

# Test
unif_test(data = hyp_data, type = "Rayleigh", p_value = "asymp")

## ----venus--------------------------------------------------------------------
# Load spherical data
data(venus)
head(venus)
nrow(venus)

# Compute Cartesian coordinates from polar coordinates
venus$X <- cbind(cos(venus$theta) * cos(venus$phi),
                 sin(venus$theta) * cos(venus$phi),
                 sin(venus$phi))

# Test uniformity using the Projected Cramér-von Mises and Projected
# Anderson-Darling tests on S^2, both using asymptotic distributions
unif_test(data = venus$X, type = c("PCvM", "PAD"), p_value = "asymp")

# Define a handler for progressr
require(progress)
require(progressr)
handlers(handler_progress(
  format = paste("(:spin) [:bar] :percent Iter: :current/:total Rate:",
                 ":tick_rate iter/sec ETA: :eta Elapsed: :elapsedfull"),
  clear = FALSE))

# Test uniformity using Monte-Carlo approximated null distributions
with_progress(
  unif_test(data = venus$X, type = c("PCvM", "PAD"),
            p_value = "MC", chunks = 1e2, M = 5e2, cores = 2)
)

## ----unif_stat_vec------------------------------------------------------------
# M samples of size n on S^2
samps_sph <- r_unif_sph(n = 30, p = 3, M = 10)

# Apply all statistics to the M samples
unif_stat(data = samps_sph, type = "all")

## ----MC-----------------------------------------------------------------------
# Break the simulation in 10 chunks of tasks to be divided between 2 cores
sim <- unif_stat_MC(n = 30, type = "all", p = 3, M = 1e2, cores = 2,
                    chunks = 10)

# Critical values for test statistics
sim$crit_val_MC

# Test statistics
head(sim$stats_MC)

# Power computation using a pre-built sampler for the alternative
pow <- unif_stat_MC(n = 30, type = "all", p = 3, M = 1e2, cores = 2,
                    chunks = 10, r_H1 = r_alt, crit_val = sim$crit_val_MC,
                    alt = "vMF", kappa = 1)
pow$power_MC

# How to use a custom sampler according to ?unif_stat_MC
r_H1 <- function(n, p, M, l = 1) {

  samp <- array(dim = c(n, p, M))
  for (j in 1:M) {

    samp[, , j] <- mvtnorm::rmvnorm(n = n, mean = c(l, rep(0, p - 1)),
                                    sigma = diag(rep(1, p)))
    samp[, , j] <- samp[, , j] / sqrt(rowSums(samp[, , j]^2))

  }
  return(samp)

}
pow <- unif_stat_MC(n = 30, type = "all", p = 3, M = 1e2, cores = 2,
                    chunks = 5, r_H1 = r_H1, crit_val = sim$crit_val_MC)
pow$power_MC

## ----crit_val-----------------------------------------------------------------
# Using precomputed critical values
ht1 <- unif_test(data = samps_sph[, , 1], type = "Rayleigh",
                 p_value = "crit_val", crit_val = sim$crit_val_MC)
ht1
ht1$reject

# Using precomputed Monte Carlo statistics
ht2 <- unif_test(data = samps_sph[, , 1], type = "Rayleigh",
                 p_value = "MC", stats_MC = sim$stats_MC)
ht2
ht2$reject

# Faster than
unif_test(data = samps_sph[, , 1], type = "Rayleigh", p_value = "MC")

## -----------------------------------------------------------------------------
citation("sphunif")

