library(RBesT)
library(parallel)
library(ggplot2)
library(latex2exp)
library(ggh4x)

## ===========================================================================================================
## RMSE computation
## ===========================================================================================================

w <- c(0.25, 0.5, 0.75)                                         ### prior weight w in paper given to historical data
n_robust <- c(0.0025, 0.04, 0.25, 1, 2)                         ## dispersion of robust component in terms of effective sample size
robust_component_location <- c("ext", "current_observed_mean") ## location of robust component
mu_ext <- seq(0, 3, by = 0.01)                                 ## prior data conflict in paper, x-axis of Figure 1

grid <- expand.grid(w = w, mu_ext = mu_ext, n_robust = n_robust, robust_component_location = robust_component_location)

true_current_mu <- 0                                             ## true mean of current data
n <- 20                                                          ## current sample size
n_ext <- 15                                                     ## historical data sample size
SIGMA <- 1                                                       # current sigma, assumed as known...
sigma <- 1                                                       # historical sigma, assumed as known...

nsims <- 1e6
ncores <- 64

MSE_results <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
  function(params) {
    with(params, {
      set.seed(123)
      y <- rnorm(nsims, mean = true_current_mu, sd = SIGMA / sqrt(n))

      if (robust_component_location == "current_observed_mean") {
        MSE <- mean(sapply(1:length(y), function(z) {
          mixprior <- mixcombine(mixnorm(informative = c(1, c(mu_ext, sigma / sqrt(n_ext))), sigma = sigma), mixnorm(vague = c(1, c(y[z], n_robust)), sigma = sigma, param = "mn"), weight = c(w, (1 - w)))
          suppressMessages({
            posterior.sum <- postmix(mixprior, m = y[z], n = n)
          })

          post_mean <- (posterior.sum["m", 1] * posterior.sum["w", 1]) + (posterior.sum["m", 2] * posterior.sum["w", 2])
          return((post_mean - true_current_mu)^2)
        }))
      } else if (robust_component_location == "ext") {
        MSE <- mean(sapply(1:length(y), function(z) {
          mixprior <- mixcombine(mixnorm(informative = c(1, c(mu_ext, sigma / sqrt(n_ext))), sigma = sigma), mixnorm(vague = c(1, c(mu_ext, n_robust)), sigma = sigma, param = "mn"), weight = c(w, (1 - w)))
          suppressMessages({
            posterior.sum <- postmix(mixprior, m = y[z], n = n)
          })

          post_mean <- (posterior.sum["m", 1] * posterior.sum["w", 1]) + (posterior.sum["m", 2] * posterior.sum["w", 2])
          return((post_mean - true_current_mu)^2)
        }))
      }
      return(data.frame(w, mu_ext, n_robust, MSE, MSE_ratio = MSE / (SIGMA / n), MSE_wo = SIGMA / n, standardized_RMSE = sqrt(MSE) / sqrt(SIGMA / n), robust_component_location))
    })
  },
  mc.cores = ncores
))

# =======================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------------#
# RMSE with a student-t robust component and the student-t is approximated by a mixture of normals (in Figure 1b in the paper)
# =======================================================================================================================#

source("student_t_prior.R")                          ## student-t prior as a mixture of normals

Nc <- c(100)                                         # number of components approximating the student-t
nu <- c(3)                                           # degrees of freedom of the student-t

grid <- expand.grid(w = w, mu_ext = mu_ext, Nc = Nc, nu = nu, n = n, n_ext = n_ext)

MSE_results_studentt <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
  function(params) {
    with(params, {
      set.seed(123)

      m0 <- mu_ext
      s0 <- 1.0 # scale of the student-t

      mix_prior_q <- mixnorm_t_prior(nu, m0, s0, nodes = Nc, method = "quantile", sigma_ref = 1)

      y <- rnorm(nsims, mean = true_current_mu, sd = SIGMA / sqrt(n))


      ## prior
      prior_mix <- tryCatch(
        {
          mixcombine(mixnorm(informative = c(1, c(mu_ext, sigma / sqrt(n_ext))), sigma = 1), mix_prior_q, weight = c(w, 1 - w))
        },
        error = function(e) {
          # Return NA if an error occurs
          NA
        }
      )
      squared_error <- sapply(1:length(y), function(z) {
        tryCatch(
          {
            (tryCatch(
              {
                suppressMessages(post_mix <- postmix(prior_mix, m = y[z], n = n))
                sum(sapply(1:ncol(post_mix), function(z) {
                  post_mix["w", z] * post_mix["m", z]
                }))
              },
              error = function(e) {
                # Return NA if an error occurs
                NA
              }
            ) - true_current_mu)^2
          },
          error = function(e) {
            # Return NA if an error occurs
            NA
          }
        )
      })
      return(data.frame(mu_ext = mu_ext, w = w, MSE = mean(squared_error, na.rm = T), MSE_ratio = mean(squared_error, na.rm = T) / (SIGMA / n), MSE_wo = SIGMA / n, standardized_RMSE = sqrt(mean(squared_error, na.rm = T)) / sqrt(SIGMA / n), t_scale = s0, Nc = Nc, df = nu, n = n, n_ext = n_ext))
    })
  },
  mc.cores = ncores
))

#-------------------------------------------------------------------------------------
## plot (Figure 1b in the paper)
#-------------------------------------------------------------------------------------
ggplot() +
  geom_hline(yintercept = 1, size = 1.2) +
  geom_line(data = MSE_results[MSE_results$w != 0.75&MSE_results$n_robust==1, ], aes(x = mu_ext, y = standardized_RMSE, linetype = as.factor(w), color = as.factor(robust_component_location)), size = 1.2) +
  geom_line(data = MSE_results_studentt[MSE_results_studentt$w != 0.75, ], aes(x = mu_ext, y = standardized_RMSE, linetype = as.factor(w), color = "darkorchid"), size = 1.2) +
  xlab(TeX("$\\bar{y}_{ext}-\\theta$")) +
  ylab("Standardized RMSE") +
  theme_bw() +
  scale_x_continuous(breaks = c(0, 1, 2, 3), labels = c("0", "1", "2", "3"), expand = c(0, 0)) +
  theme(
    axis.text = element_text(color = c("black", "black", "black", "black", "black")),
    axis.ticks = element_line(color = c("black", "black", "black", "black", "black"))
  ) +
  theme(
    axis.text = element_text(size = 28, family = "serif"),
    axis.title.x = element_text(size = 28, family = "serif"),
    axis.title.y = element_text(size = 30, family = "serif"),
    legend.text = element_text(size = 28, family = "serif"),
    legend.title = element_text(size = 30, family = "serif"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.width = unit(6, "line"),
    legend.spacing.x = unit(1.2, "cm"),
    legend.spacing.y = unit(-0.1, "cm")
  ) +
  guides(color = "none", linetype = "none") +
  scale_color_manual(values = c("darkgreen", "red", "darkorchid"), labels = unname(TeX(c("$\\mu_{robust}=\\bar{y}$", "$\\mu_{robust}=\\bar{y}_{ext}$", "$\\mu_{robust}=\\bar{y}_{ext} \\, (t \\, distr.)$")))) +
  theme(panel.background = element_rect(fill = "white", colour = NA)) +
  theme(
    axis.line = element_line(colour = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )
