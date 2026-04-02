library(RBesT)
library(parallel)

## ===========================================================================================================
## TIE rate and power in a restricted drift range (Table 1 in the paper)
## ===========================================================================================================

source("power.normal.R")
source("TIErate_hybridcontrol.R")
source("power_hybridcontrol.R")

w <- c(0.5)                       # weight given to the informative component (0 < weight < 1).
mu_ext <- 0                      ##  true mean in historical control
n_ext <- 15                      # external control sample size
n_c <- c(20)                      # current control sample size
n_t <- c(20)                      # current treatment arm sample size
n_robust <- c(1)                  # number of observations the non-informative prior corresponds to
robust_component_location <- c("ext", "current_observed_mean") ## location of robust component
delta_index <- 1:4
grid <- expand.grid(w = w, mu_ext = mu_ext, robust_component_location = robust_component_location, n_ext = n_ext, n_c = n_c, n_t = n_t, n_robust = n_robust, delta_index = delta_index)

TREAT.SIGMA <- 1                  # current sigma in treatment arm(assumed known)
CONTROL.SIGMA <- 1                # current control sigma (assumed known)
H1Diff <- 0.8330874
cE <- .975

nsims <- 1e6
ncores <- 64

restricted_drift_OCs <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
  function(params) {
    with(params, {
      H1Diff <- 0.8330874                      # where to evaluate power: theta_c - theta_t=0.8330874
      wrapper_func_TIErate <- function(vec) {
        TIErate_hybridcontrol(
          H0Diff = 0,                          # theta_c - theta_t=0
          TREAT.SIGMA = 1,                     # current sigma in treatment arm(assumed known)
          CONTROL.SIGMA = 1,                   # current control sigma (assumed known)
          n_ext = n_ext,                       # external control sample size
          n_c = n_c,                           # current control sample size
          n_t = n_t,                           # current treatment sample size
          ext.SIGMA = 1,                       # external control sigma (assumed known)
          n_robust = n_robust,                 # number of observations the non-informative prior corresponds to
          muE = 0, cE = .975,                  # P(treat-control>muE|.)> cE
          mu_ext = mu_ext,                     ##  true mean in historical control
          Delta = vec,                         # this is difference between observed mean historical control and theta_c divided by sigma (i.e. value on x-axis)
          w = w,                               # weight given to the informative component (0 < weight < 1).
          robust_component_location = robust_component_location, # mean of the non-informative component.
          nsim = nsims                         ## number of monte carlo simulations
        )
      }

      wrapper_func_power <- function(vec) {
        power_hybridcontrol(
          H1Diff = 0.8330874,                 # where to evaluate power: theta_c - theta_t=0.8330874
          TREAT.SIGMA = 1,                     # current sigma in treatment arm(assumed known)
          CONTROL.SIGMA = 1,                   # current control sigma (assumed known)
          n_ext = n_ext,                       # external control sample size
          n_c = n_c,                           # current control sample size
          n_t = n_t,                           # current treatment sample size
          ext.SIGMA = 1,                       # external control sigma (assumed known)
          n_robust = n_robust,                 # number of observations the non-informative prior corresponds to
          muE = 0, cE = .975,                  # P(treat-control>muE|.)> cE
          mu_ext = mu_ext,                     ##  true mean in historical control
          Delta = vec,                         # this is difference between observed mean historical control and theta_c divided by sigma (i.e. value on x-axis)
          w = w,                               # weight given to the informative component (0 < weight < 1).
          robust_component_location = robust_component_location, # mean of the non-informative component.
          nsim = nsims                        ## number of monte carlo simulations
        )
      }

      delta_lower <- c(-0.1, -0.2, -0.4, -0.5)[delta_index]          ## drift restricted by delta
      delta_upper <- c(0.1, 0.2, 0.4, 0.5)[delta_index]              ## drift restricted by delta

      max_TIE_in_interval <- optimize(wrapper_func_TIErate, interval = c(delta_lower, delta_upper), maximum = T)$objective

      powerwo <- Powerwo_2sample(treat.mu = H1Diff, control.mu = 0, treat.sigma = TREAT.SIGMA, control.sigma = CONTROL.SIGMA, n = n_t, alpha = 1 - cE)

      # power w/o borrowing callibrated to borrowing:
      powerwo_alphaBMix <- Powerwo_2sample(treat.mu = H1Diff, control.mu = 0, treat.sigma = TREAT.SIGMA, control.sigma = CONTROL.SIGMA, n = n_t, alpha = max_TIE_in_interval)

      max_powerwMix <- optimize(wrapper_func_power, interval = c(delta_lower, delta_upper), maximum = T)$objective

      max_powergain <- max_powerwMix - powerwo_alphaBMix

      return(data.frame(mu_ext = mu_ext, w = w, n_robust, robust_component_location, n_ext, n_t, n_c, abs_delta = delta_upper, max_TIE_in_interval, powerwo, powerwo_alphaBMix, max_powerwMix, max_powergain))
    })
  },
  mc.cores = ncores
))
