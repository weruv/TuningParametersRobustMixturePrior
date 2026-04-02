library(RBesT)
library(parallel)
library(ggplot2)
library(latex2exp)
library(ggh4x)

## ===========================================================================================================
## TIE rate and power- one arm
## ===========================================================================================================

source("function_1arm.R")

w <- c(0.25, 0.5, 0.75)                                          ### prior weight w in paper given to historical data
n_robust <- c(0.0025, 0.04, 0.25, 1, 2)                          ## dispersion of robust component in terms of effective sample size
robust_comp_location <- c("ext", "H0", "current_data_mean")     ## location of robust component
mu_ext <- seq(-3, 3, by = 0.01)                                 ## prior data conflict in paper, x-axis of Figure 1

grid <- expand.grid(w = w, n_robust = n_robust, robust_comp_location = robust_comp_location, mu_ext = mu_ext)

nsims <- 1e6
ncores <- 64

TIErate_n_power_1arm <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
  function(params) {
    with(params, {
      RhpcBLASctl::blas_set_num_threads(1)


      alpha <- 0.025

      theta1 <- 0.5            # where I evaluate power, i.e. mean for current
      H0MU <- 0                # where I evaluate type I error rate

      N <- 20                  # sample size current data (n in paper)
      SIGMA <- 1               # current sigma, assumed as known...

      muE <- 0                 #    P(x>muE|.)> cE
      cE <- .975
      sigma <- 1               # historical sigma, assumed as known...
      n_ext <- 15              # sample size historical data (n_ext in paper)


      results_inc <- rep(NA, 6)

      for (dE in mu_ext) # grid of mean mu_ext

      {
        info <- c(dE, sigma / sqrt(n_ext))

        set.seed(123)
        y <- rnorm(nsims, mean = H0MU, sd = SIGMA / sqrt(N))

        decision <- sapply(1:length(y), function(z) {
          pH1.onearm(
            yt = y[z], nt = N, priorPars = info, w.mix = w, sigma = 1,
            bound = FALSE, pH1array = NULL, a0t = 0.5, a0c = 0.5, n0 = 1, th0 = muE, sigma_robust = 1 / sqrt(n_robust), robust_location = robust_comp_location
          )
        })
        alphaB <- mean(decision > cE)

        set.seed(123)
        y <- rnorm(nsims, mean = theta1, sd = SIGMA / sqrt(N))

        decision_h1 <- sapply(1:length(y), function(z) {
          pH1.onearm(
            yt = y[z], nt = N, priorPars = info, w.mix = w, sigma = 1,
            bound = FALSE, pH1array = NULL, a0t = 0.5, a0c = 0.5, n0 = 1, th0 = muE, sigma_robust = 1 / sqrt(n_robust), robust_location = robust_comp_location
          )
        })
        powerwMix <- mean(decision_h1 > cE)

        # #compare power for w and w/o borrowing, taking the adjusted typeI error
        powerwo_alphaB <- pnorm(sqrt(N) * (theta1 - H0MU) / SIGMA - qnorm(1 - alphaB)) # frequentist power for adjusted alpha
        powerwo_unadjusted <- pnorm(sqrt(N) * (theta1 - H0MU) / SIGMA - qnorm(cE)) # frequentist power for unadjusted alpha

        powerwMixdiff <- powerwMix - powerwo_alphaB
        results_inc <- rbind(results_inc, as.numeric(c(dE, alphaB, powerwo_alphaB, powerwMix, powerwMixdiff, powerwo_unadjusted)))
      }

      resultsrob <- data.frame(results_inc)
      resultsrob <- resultsrob[-1, ]
      names(resultsrob) <- c("mu_ext", "alphaB", "powerwo", "powerwMix", "powerwMixdiff", "powerwo_unadjusted")
      return(cbind(resultsrob, prior_weight = w, n_curr = N, n_ext = n_ext, n_robust = n_robust, robust_comp_location = robust_comp_location))
    })
  },
  mc.cores = ncores
))

# =======================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------------#
# TIE rate with a student-t robust component and the student-t is approximated by a mixture of normals (in Figure 1a in the paper)
# =======================================================================================================================#

source("student_t_prior.R")                         ## student-t prior as a mixture of normals

N <- 20                                             # sample size current data (n in paper)
n_ext <- c(15)

Nc <- c(100)                                        # number of components approximating the student-t
nu <- c(3)                                          # degrees of freedom of the student-t

grid <- expand.grid(w = w, mu_ext = mu_ext, Nc = Nc, nu = nu, N = N, n_ext = n_ext)

TIErate_onearm_studentt <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
  function(params) {
    with(params, {
      set.seed(123)


      H0MU <- 0             # where I evaluate type I error rate

      SIGMA <- 1            # current sigma, assumed as known...
      sigma <- 1            # historical sigma, assumed as known...

      m0 <- mu_ext
      s0 <- 1.0             # scale of the student-t

      mix_prior_q <- mixnorm_t_prior(nu, m0, s0, nodes = Nc, method = "quantile", sigma_ref = 1)

      y <- rnorm(nsims, mean = H0MU, sd = SIGMA / sqrt(N))

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
      decision <- sapply(1:length(y), function(z) {
        post.mix <- tryCatch(
          {
            suppressMessages(postmix(prior_mix, m = y[z], n = N))
          },
          error = function(e) {
            # Return NA if an error occurs
            NA
          }
        )
        post_p <- tryCatch(
          {
            pmix(
              post.mix,
              q = H0MU,
              lower.tail = TRUE
            ) <= 0.025
          },
          error = function(e) {
            # Return NA if an error occurs
            NA
          }
        )
        return(post_p)
      })
      return(data.frame(mu_ext = mu_ext, w = w, TIE_rate = mean(decision, na.rm = T), t_scale = s0, Nc = Nc, df = nu, N = N, n_ext = n_ext))
    })
  },
  mc.cores = ncores
))

#-------------------------------------------------------------------------------------
## plot (Figure 1a in the paper)
#-------------------------------------------------------------------------------------
ggplot() +
  geom_hline(yintercept = 0.025, size = 1.2) +
  geom_line(data = TIErate_n_power_1arm[TIErate_n_power_1arm$prior_weight != 0.75 & TIErate_n_power_1arm$n_robust == 1, ], aes(x = mu_ext, y = alphaB, linetype = as.factor(prior_weight), color = as.factor(robust_comp_location)), size = 1.2) +
  geom_line(data = TIErate_onearm_studentt[TIErate_onearm_studentt$w != 0.75, ], aes(x = mu_ext, y = TIE_rate, linetype = as.factor(w), color = "darkorchid"), size = 1.2) +
  xlab(TeX("$\\bar{y}_{ext}-\\theta_0$")) +
  ylab("Type I error rate") +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.1), labels = c("0", "0.025", "0.05", "0.10"), limits = c(0, 0.10)) +
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3), labels = c("-3", "-2", "-1", "0", "1", "2", "3"), limits = c(-3, 3), expand = c(0, 0)) +
  theme(
    axis.text = element_text(color = c("black", "black", "black", "black", "black")),
    axis.ticks = element_line(color = c("black", "black", "black", "black", "black"))
  ) +
  theme(
    axis.text = element_text(size = 28, family = "serif"),
    axis.title.x = element_text(size = 28, family = "serif"),
    axis.title.y = element_text(size = 30, family = "serif"),
    legend.text = element_text(size = 30, family = "serif"),
    legend.title = element_text(size = 28, family = "serif"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.width = unit(6, "line"),
    legend.spacing.x = unit(1.2, "cm"),
    legend.spacing.y = unit(-0.1, "cm")
  ) +
  guides(color = "none", linetype = "none") +
  scale_color_manual(values = c("darkgreen", "deepskyblue", "red", "darkorchid"), labels = unname(TeX(c("$\\mu_{robust}=\\bar{y}$", "$\\mu_{robust}=\\theta_0$", "$\\mu_{robust}=\\bar{y}_{ext}$", "$\\mu_{robust}=\\bar{y}_{ext} \\, (t \\, distr.)$")))) +
  theme(panel.background = element_rect(fill = "white", colour = NA)) +
  theme(
    axis.line = element_line(colour = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )