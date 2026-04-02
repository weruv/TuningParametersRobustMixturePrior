library(RBesT)
library(parallel)
library(ggplot2)
library(latex2exp)
library(ggh4x)

## ==================================================================================================
## TIE rate and power for the hybrid control design 
## ==================================================================================================
source("power.normal.R")
source("TIErate_n_power_hybridcontrol.R")

w <- c(0.25, 0.5, 0.75)                              # weight given to the informative component (0 < weight < 1).
mu_ext <- 0                                         ##  true mean in historical control
Delta <- seq(-2, 2, by = 0.01)                       # this is difference between observed mean historical control and theta_c divided by sigma (i.e. value on x-axis)
n_ext <- 15                                         # external control sample size
n_c <- 20                                            # current control sample size
n_t <- 20                                            # current treatment arm sample size
n_robust <- 1                                        # number of observations the non-informative prior corresponds to
robust_component_location <- c("ext", "current_observed_mean") ## location of robust component

grid <- expand.grid(w = w, mu_ext = mu_ext, Delta = Delta, robust_component_location = robust_component_location, n_ext = n_ext, n_c = n_c, n_t = n_t, n_robust = n_robust)
ncores <- 64
nsims <- 1e6

TIErate_n_power_hybridcontrol <- do.call(rbind, mclapply(1:nrow(grid), function(index) {
  TIErate_n_power_hybridcontrol(
    H0Diff = 0,                                       # theta_c - theta_t=0
    H1Diff = 0.8330874,                               # where to evaluate power: theta_c - theta_t=0.8330874
    TREAT.SIGMA = 1,                                  # current sigma in treatment arm(assumed known)
    CONTROL.SIGMA = 1,                                # current control sigma (assumed known)
    n_ext = grid$n_ext[index],                      # external control sample size
    n_c = grid$n_c[index],                            # current control sample size
    n_t = grid$n_t[index],                            # current treatment sample size
    ext.SIGMA = 1,                                   # external control sigma (assumed known)
    n_robust = grid$n_robust[index],                  # number of observations the non-informative prior corresponds to
    muE = 0, cE = .975,                               # P(treat-control>muE|.)> cE
    mu_ext = grid$mu_ext[index],                    ##  true mean in historical control
    Delta = grid$Delta[index],                        # this is difference between observed mean historical control and theta_c divided by sigma (i.e. value on x-axis)
    w = grid$w[index],                                # weight given to the informative component (0 < weight < 1).
    robust_component_location = grid$robust_component_location[index], # mean of the non-informative component.
    nsim = nsims                                      ## number of monte carlo simulations
  )
}, mc.cores = ncores))

## ------------------------------------------------------------------------------------------------
## Max TIE rate 
## the code below computes the max TIE rate which is then used to calibrate the power without borrowing
## ------------------------------------------------------------------------------------------------
source("TIErate_hybridcontrol.R")

robust_component_location <- c("current_observed_mean") ## location of robust component, the max TIE for robust component at location "ext" is 1

grid <- expand.grid(w = w, mu_ext = mu_ext, robust_component_location = robust_component_location, n_ext = n_ext, n_c = n_c, n_t = n_t, n_robust = n_robust)

max_TIErate <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
  function(params) {
    with(params, {
      wrapper_func <- function(vec) {
        TIErate_hybridcontrol(
          H0Diff = 0,                   # theta_c - theta_t=0
          TREAT.SIGMA = 1,              # current sigma in treatment arm(assumed known)
          CONTROL.SIGMA = 1,            # current control sigma (assumed known)
          n_ext = n_ext,                # external control sample size
          n_c = n_c,                    # current control sample size
          n_t = n_t,                    # current treatment sample size
          ext.SIGMA = 1,                # external control sigma (assumed known)
          n_robust = n_robust,          # number of observations the non-informative prior corresponds to
          muE = 0, cE = .975,           # P(treat-control>muE|.)> cE
          mu_ext = mu_ext,              ##  true mean in historical control
          Delta = vec,                  # this is difference between observed mean historical control and theta_c divided by sigma (i.e. value on x-axis)
          w = w,                        # weight given to the informative component (0 < weight < 1).
          robust_component_location = robust_component_location, # mean of the non-informative component.
          nsim = nsims                  ## number of monte carlo simulations
        )
      }
      max_alphaB_2arm <- optimize(wrapper_func, interval = c(-0.5, 1), maximum = T)$objective
      return(data.frame(mu_ext = mu_ext, w = w, n_robust, robust_component_location, n_ext, n_t, n_c, max_alphaB_2arm))
    })
  },
  mc.cores = ncores
))

## ------------------------------------------------------------------------------------------------
## calibrate power without borrowing to max TIE rate
## ------------------------------------------------------------------------------------------------
max_TIErate$power_calibrated <- Powerwo_2sample(treat.mu = 0.8330874, control.mu = 0, treat.sigma = 1, control.sigma = 1, n = n_c, alpha = max_TIErate$max_alphaB_2arm)

## ------------------------------------------------------------------------------------------------
## merge power_calibrated to rest of dataset
TIErate_n_power_hybridcontrol <- merge(TIErate_n_power_hybridcontrol, max_TIErate, by = colnames(max_TIErate)[!colnames(max_TIErate) %in% c("max_alphaB_2arm", "power_calibrated")], all.x = T)

## for robust component with unit-information located at ext, max TIE rate is 1
TIErate_n_power_hybridcontrol$max_alphaB_2arm[TIErate_n_power_hybridcontrol$robust_component_location == "ext"] <- 1
TIErate_n_power_hybridcontrol$power_calibrated[TIErate_n_power_hybridcontrol$robust_component_location == "ext"] <- 1

## compute the difference between power with borrowing and power calibrated to borrowing (for Figure A.11)
TIErate_n_power_hybridcontrol$powermixdiff <- TIErate_n_power_hybridcontrol$powerWMix - TIErate_n_power_hybridcontrol$power_calibrated

### plot (Figure 5 in the paper)
TIErate_n_power_hybridcontrol$robust_component_location <- as.factor(TIErate_n_power_hybridcontrol$robust_component_location)
levels(TIErate_n_power_hybridcontrol$robust_component_location) <- c(current_observed_mean = latex2exp::TeX("$\\mu_{robust}=\\bar{y}_{c}$"), ext = latex2exp::TeX("$\\mu_{robust}=\\bar{y}_{ext}$"))
TIErate_n_power_hybridcontrol$title <- "Location"

ggplot() +
  geom_hline(yintercept = c(0.025, 0.75), size = 1.1, linetype = "dotdash") +
  geom_line(data = TIErate_n_power_hybridcontrol[TIErate_n_power_hybridcontrol$w != .75, ], aes(x = Delta, y = powerWMix, color = as.factor(w), linetype = "Power with borrowing"), size = 1.1) +
  geom_line(data = TIErate_n_power_hybridcontrol[TIErate_n_power_hybridcontrol$w != .75, ], aes(x = Delta, y = TIE_rate, color = as.factor(w), linetype = "TIE rate"), size = 1.1) +
  geom_line(data = TIErate_n_power_hybridcontrol[TIErate_n_power_hybridcontrol$w != .75, ], aes(x = Delta, y = power_calibrated, color = as.factor(w), linetype = "Power calibrated \n to borrowing"), size = 1.1) +
  xlab(TeX("$\\theta_c-\\bar{y}_{ext}$")) +
  ylab(TeX("$\\Probability \\,to \\,reject \\,H_0")) +
  scale_linetype_manual(values = c("Power with borrowing" = "solid", "TIE rate" = "longdash", "Power calibrated \n to borrowing" = "dotted")) +
  scale_y_continuous(breaks = c(0, 0.025, 0.06, 0.60, 0.75, 0.90, 1), labels = c("0", "0.025", "0.06", "0.60", "0.75", "0.90", "1")) +
  scale_x_continuous(breaks = c(-1, 0, 1), labels = c("-1", "0", "1"), expand = c(0, 0)) +
  ggbreak::scale_y_break(c(0.06, 0.58), scales = 1.5, space = 0.2) +
  facet_wrap(. ~ title * robust_component_location, labeller = labeller(.cols = label_parsed), strip.position = "top", scales = "free_x") + ### this makes the two plots next to each other
  theme_bw() +
  theme(
    axis.text = element_text(color = c("black", "black", "black", "black", "black")),
    axis.ticks = element_line(color = c("black", "black", "black", "black", "black"))
  ) +
  theme(
    axis.text.x = element_text(size = 28, family = "serif"),
    axis.text.y = element_text(size = 28, family = "serif"),
    axis.title.x = element_text(size = 28, family = "serif"),
    axis.title.y = element_text(size = 28, family = "serif"),
    legend.text = element_text(size = 26, family = "serif"),
    legend.title = element_text(size = 26, family = "serif"),
    legend.position = "right",
    legend.key.width = unit(3.5, "line"),
    strip.text.x = element_text(size = 28, family = "serif"),
    strip.text.y = element_text(size = 28, family = "serif")
  ) +
  labs(color = "w", linetype = "") +
  theme(panel.background = element_rect(fill = "white", colour = NA)) +
  theme(
    axis.line = element_line(colour = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  ) +
  theme(
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.line.x.top = element_blank(),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.line.y.right = element_blank(),
    legend.position = "right",
    legend.margin = margin(10, 0, 0, 0),
    legend.box.margin = margin(t = 10)
  )

#=======================================================================================================================#
#-----------------------------------------------------------------------------------------------------------------------#
# TIE rate with a student-t robust component and the student-t is approximated by a mixture of normals (Figure 6 in the paper) 
#=======================================================================================================================#

source("student_t_prior.R")                         ## student-t prior as a mixture of normals

mu_ext <- seq(-2, 2, by = 0.01)                        # mean historical control

Nc <- c(100)                                       # number of components approximating the student-t
nu <- c(3)                                         # degrees of freedom of the student-t

grid <- expand.grid(w = w, mu_ext = mu_ext, Nc = Nc, nu = nu, n_c = n_c, n_t = n_t, n_ext = n_ext)

TIErate_hybridcontrol_studentt <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
  function(params) {
    with(params, {
      set.seed(123)
      
      TREAT.SIGMA <- 1             # current sigma in treatment arm(assumed known)
      CONTROL.SIGMA <- 1           # current control sigma (assumed known)
      ext.SIGMA <- 1               # external control sigma (assumed known)

      H0Diff <- 0

      m0 <- mu_ext
      s0 <- 1.0                    # scale of the student-t

      mix_prior_q <- mixnorm_t_prior(nu, m0, s0, nodes = Nc, method = "quantile", sigma_ref = 1)

      dat_control <- rnorm(nsims, mean = H0Diff, sd = CONTROL.SIGMA / sqrt(n_c))
      dat_treat <- rnorm(nsims, mean = H0Diff, sd = TREAT.SIGMA / sqrt(n_t))

      ## prior for the control arm
      prior_control <- tryCatch(
        {
          mixcombine(mixnorm(informative = c(1, c(mu_ext, ext.SIGMA / sqrt(n_ext))), sigma = 1), mix_prior_q, weight = c(w, 1 - w))
        },
        error = function(e) {
          # Return NA if an error occurs
          NA
        }
      )

      decision <- sapply(1:length(dat_control), function(z) {
        ## posterior for the control arm
        post.control <- tryCatch(
          {
            suppressMessages(postmix(prior_control, m = dat_control[z], n = n_c))
          },
          error = function(e) {
            # Return NA if an error occurs
            NA
          }
        )
        ## posterior for the treatment arm
        post.treat <- mixnorm(
          comp1 = c(1, dat_treat[z], TREAT.SIGMA / sqrt(n_t)),
          sigma = 1
        ) # assumes flat prior
        post_p <- tryCatch(
          {
            pmixdiff(
              mix1 = post.treat,
              mix2 = post.control,
              q = 0,
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
      return(data.frame(dE = mu_ext, w = w, TIE_rate = mean(decision, na.rm = T), t_scale = s0, Nc = Nc, df = nu, n_c = n_c, n_t = n_t, n_ext = n_ext))
    })
  },
  mc.cores = ncores
))

## plot (Figure 6) ---------------------------------------------------------------------------------------------
ggplot() +
  geom_hline(aes(yintercept = c(0.025)), size = 1.1) +
  geom_line(data = TIErate_n_power_hybridcontrol[TIErate_n_power_hybridcontrol$w == .5 & TIErate_n_power_hybridcontrol$robust_component_location == "ext", ], aes(x = Delta, y = TIE_rate, color = "red"), size = 1.1) +
  geom_line(data = TIErate_hybridcontrol_studentt[TIErate_hybridcontrol_studentt$w == 0.5, ], aes(x = 0 - dE, TIE_rate, color = "darkorchid"), size = 1.1) +
  xlab(TeX("$\\theta_c-\\bar{y}_{ext}$")) +
  ylab("Type I error rate") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    axis.text = element_text(color = c("black", "black", "black", "black", "black")),
    axis.ticks = element_line(color = c("black", "black", "black", "black", "black"))
  ) +
  theme(
    axis.text = element_text(size = 28, family = "serif"),
    axis.title.x = element_text(size = 28, family = "serif"),
    axis.title.y = element_text(size = 28, family = "serif"),
    legend.text = element_text(size = 26, family = "serif"),
    legend.title = element_text(size = 26, family = "serif"),
    legend.position = "bottom",
    legend.text.align = 0,
    legend.box = "vertical",
    legend.key.width = unit(6, "line"),
    legend.spacing.x = unit(1.2, "cm"),
    legend.spacing.y = unit(-0.1, "cm")
  ) +
  guides(color = guide_legend(title = "")) +
  scale_color_manual(values = c("darkorchid" = "darkorchid", "red" = "red"), labels = c("red" = TeX("$\\mu_{robust}=\\bar{y}_{ext}  \\, (normal \\, distr.)$"), "darkorchid" = TeX("$\\mu_{robust}=\\bar{y}_{ext} \\, (t \\, distr.)$"))) +
  theme(panel.background = element_rect(fill = "white", colour = NA)) +
  theme(
    axis.line = element_line(colour = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )