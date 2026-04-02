library(RBesT)
library(parallel)
library(RColorBrewer)
library(tidyr)
library(latex2exp)
library(ggh4x)

## ===========================================================================================================
## sweet spot computation
## ===========================================================================================================

source("power.normal.R")
source("TIErate_hybridcontrol.R")
source("power_hybridcontrol.R")

w <- seq(0, 1, by = 0.07)                              # weight given to the informative component (0 < weight < 1).
mu_ext <- 0                                           ##  true mean in historical control
n_ext <- 15                                           # external control sample size
n_c <- c(10, 20, 30)                                   # current control sample size
n_t <- c(10, 20, 30)                                   # current treatment arm sample size
n_robust <- c(0.0025, 0.04, 0.25, 1)                   # number of observations the robust component corresponds to
robust_component_location <- c("ext", "current_observed_mean") ## location of robust component

grid <- expand.grid(w = w, mu_ext = mu_ext, robust_component_location = robust_component_location, n_ext = n_ext, n_c = n_c, n_t = n_t, n_robust = n_robust)

TREAT.SIGMA <- 1                                       # current sigma in treatment arm(assumed known)
CONTROL.SIGMA <- 1                                     # current control sigma (assumed known)
H1Diff <- 0.8330874
cE <- .975

nsims <- 1e6
ncores <- 64

sweetspot <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
  function(params) {
    with(params, {
      wrapper_func_TIErate <- function(vec) {
        TIErate_hybridcontrol(
          H0Diff = 0,                        # theta_c - theta_t=0
          TREAT.SIGMA = 1,                   # current sigma in treatment arm(assumed known)
          CONTROL.SIGMA = 1,                 # current control sigma (assumed known)
          n_ext = n_ext,                     # external control sample size
          n_c = n_c,                         # current control sample size
          n_t = n_t,                         # current treatment sample size
          ext.SIGMA = 1,                     # external control sigma (assumed known)
          n_robust = n_robust,               # number of observations the non-informative prior corresponds to
          muE = 0, cE = .975,                # P(treat-control>muE|.)> cE
          mu_ext = mu_ext,                   ##  true mean in historical control
          Delta = vec,                       # this is difference between observed mean historical control and theta_c divided by sigma (i.e. value on x-axis)
          w = w,                             # weight given to the informative component (0 < weight < 1).
          robust_component_location = robust_component_location, # mean of the non-informative component.
          nsim = nsims                       ## number of monte carlo simulations
        )
      }
      wrapper_func_power <- function(vec) {
        power_hybridcontrol(
          H1Diff = 0.8330874,              # where to evaluate power: theta_c - theta_t=0.8330874
          TREAT.SIGMA = 1,                  # current sigma in treatment arm(assumed known)
          CONTROL.SIGMA = 1,                # current control sigma (assumed known)
          n_ext = n_ext,                    # external control sample size
          n_c = n_c,                        # current control sample size
          n_t = n_t,                        # current treatment sample size
          ext.SIGMA = 1,                    # external control sigma (assumed known)
          n_robust = n_robust,              # number of observations the non-informative prior corresponds to
          muE = 0, cE = .975,               # P(treat-control>muE|.)> cE
          mu_ext = mu_ext,                  ##  true mean in historical control
          Delta = vec,                      # this is difference between observed mean historical control and theta_c divided by sigma (i.e. value on x-axis)
          w = w,                            # weight given to the informative component (0 < weight < 1).
          robust_component_location = robust_component_location, # mean of the non-informative component.
          nsim = nsims                      ## number of monte carlo simulations
        )
      }

      powerwo <- Powerwo_2sample(treat.mu = H1Diff, control.mu = 0, treat.sigma = TREAT.SIGMA, control.sigma = CONTROL.SIGMA, n = n_t, alpha = 1 - cE)

      upper_bound <- tryCatch(stats::uniroot(\(vec) wrapper_func_TIErate(vec) - 0.025, c(-0.5, 0.7), extendInt = "yes")$root, error = function(err) NA)
      lower_bound <- tryCatch(stats::uniroot(\(vec) wrapper_func_power(vec) - powerwo, c(-0.7, 0.5), extendInt = "yes")$root, error = function(err) NA)

      max_power_sweetspot <- optimize(wrapper_func_power, interval = c(lower_bound, upper_bound), maximum = T)

      return(data.frame(mu_ext = mu_ext, w = w, n_robust, robust_component_location, n_ext, n_t, n_c, lower_bound, upper_bound, max_power_in_sweetspot = max_power_sweetspot$objective, max_power_spot = max_power_sweetspot$maximum))
    })
  },
  mc.cores = ncores
))

## ---------------------------------------------------------------------------------------------------------------------
## plot (Figure 7 in the paper)
## ---------------------------------------------------------------------------------------------------------------------
myPalette <- colorRampPalette((brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), limits = c(0.75, max(sweetspot$max_power_in_sweetspot)))

sweetspot$title <- as.factor("Robust_precision")
levels(sweetspot$title) <- c(Robust_precision = latex2exp::TeX("$\\n_{robust}\\, (precision \\, increasing \\, from \\, left \\, to \\, right)$"))
sweetspot$robust_component_location <- as.factor(sweetspot$robust_component_location)
levels(sweetspot$robust_component_location) <- c(current_observed_mean = TeX("$\\mu_{robust}=\\bar{y}_c$"), ext = TeX("$\\mu_{robust}=\\bar{y}_{ext}$"))

sweetspot %>%
  filter(!(w %in% c(0)), n_t %in% c(20), n_c %in% c(20)) %>%
  ggplot(aes(x = w)) +
  geom_linerange(aes(ymin = lower_bound, ymax = upper_bound, x = w, color = max_power_in_sweetspot),
    size = 1.5, alpha = 1
  ) +
  geom_point(aes(y = max_power_spot), size = 1.3) +
  coord_flip() +
  facet_nested("Location" * robust_component_location ~ title * n_robust, labeller = labeller(.rows = label_parsed, .cols = label_parsed)) +
  ylab(TeX("$\\theta_c-\\bar{y}_{ext}$")) +
  xlab("Prior weight") +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 1, by = 0.21)) +
  scale_y_continuous(breaks = c(-0.4, -0.2, 0.0)) +
  theme(
    axis.text = element_text(color = c("black", "black", "black", "black", "black")),
    axis.ticks = element_line(color = c("black", "black", "black", "black", "black"))
  ) +
  theme(
    axis.text = element_text(size = 28, family = "serif"),
    axis.title.x = element_text(size = 28, family = "serif"),
    axis.title.y = element_text(size = 28, family = "serif"),
    legend.text = element_text(size = 30, family = "serif"),
    legend.title = element_text(size = 30, family = "serif"),
    legend.key.height = unit(3, "line"),
    strip.text.x = element_text(size = 30, family = "serif"),
    strip.text.y = element_text(size = 30, family = "serif"),
    panel.spacing.y = unit(1.5, "lines")
  ) +
  labs(color = "Maximum \n power") +
  sc
