library(ggplot2)
library(latex2exp)
library(ggh4x)
library(RBesT)
library(parallel)

##==========================================================================================================
# Figure 3 
##==========================================================================================================
n_ext <- 15                     ## historical data sample size
n_curr <- 20                     ## current sample size

mu_ext <- c(0, 2 / sqrt(n_ext), 8 / sqrt(n_ext))         ## no, moderate, extreme drift
robust_component_sd <- c(sqrt(0.5), 1, 2, 5, 20)            ## Dispersion of robust component
w <- seq(0, 1, by = 0.01)                                   ## prior weight w in paper given to historical data

grid <- expand.grid(w = w, mu_ext = mu_ext, robust_component_sd = robust_component_sd, n_ext = n_ext, n_curr = n_curr)

SIGMA <- 1                           # current sigma, assumed as known...
sigma <- 1                           # historical sigma, assumed as known...
current.mu <- 0                      ## true mean of the current data

ncores <- 64

posterior_weight <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
  function(params) {
    with(params, {
      RhpcBLASctl::blas_set_num_threads(1)
      prior <- mixnorm(informative = c(w, c(mu_ext, sigma / sqrt(n_ext))), vague = c(1 - w, c(mu_ext, robust_component_sd)), sigma = sigma)
      posterior_weight <- suppressMessages(postmix(prior, m = current.mu, n = n_curr)["w", 1])
      return(data.frame(mu_ext, prior_weight = w, n_curr, n_ext, robust_component_sd, n_robust = 1 / robust_component_sd^2, posterior_weight))
    })
  },
  mc.cores = ncores
))

drift.labs <- c("No drift", "Moderate drift", "Extreme drift")
names(drift.labs) <- unique(posterior_weight$mu_ext)

posterior_weight$title <- "Robust_precision"
posterior_weight$title <- as.factor(posterior_weight$title)
levels(posterior_weight$title) <- c(Robust_precision = TeX("$\\n_{robust}\\, (precision \\, increasing \\, from \\, left \\, to \\, right)$"))

ggplot(data = posterior_weight, aes(x = prior_weight, y = posterior_weight)) +
  geom_point() +
  ylab("Posterior weight on informative component") +
  xlab("Prior weight on informative component") +
  facet_nested("drift" * mu_ext ~ title * n_robust, labeller = labeller(mu_ext = drift.labs, .cols = label_parsed)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(color = c("black", "black", "black", "black", "black")),
    axis.ticks.y = element_line(color = c("black", "black", "black", "black", "black"))
  ) +
  theme(
    axis.text = element_text(size = 28, family = "serif"),
    axis.title.x = element_text(size = 32, family = "serif"),
    axis.title.y = element_text(size = 32, family = "serif"),
    legend.text = element_text(size = 30, family = "serif"),
    legend.title = element_text(size = 30, family = "serif"),
    strip.text.x = element_text(size = 32, family = "serif"),
    strip.text.y = element_text(size = 32, family = "serif")
  ) +
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

##==========================================================================================================
## Figure 4
##==========================================================================================================
mu_ext <- seq(-2, 2, by = 0.01)                          ## prior data drift in paper, x-axis of Figure 1
robust_component_sd <- c(sqrt(0.5), 1, 2, 5, 20)          ## Dispersion of robust component
w <- 0.5                                                  ## prior weight w in paper given to historical data
grid <- expand.grid(w = w, mu_ext = mu_ext, robust_component_sd = robust_component_sd, n_ext = n_ext, n_curr = n_curr)

SIGMA <- 1                     # current sigma, assumed as known...
sigma <- 1                     # historical sigma, assumed as known...
current.mu <- 0                ## true mean of the current data

posterior_weight <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
  function(params) {
    with(params, {
      RhpcBLASctl::blas_set_num_threads(1)
      prior <- mixnorm(informative = c(w, c(mu_ext, sigma / sqrt(n_ext))), vague = c(1 - w, c(mu_ext, robust_component_sd)), sigma = sigma)
      posterior_weight <- suppressMessages(postmix(prior, m = current.mu, n = n_curr)["w", 1])
      return(data.frame(mu_ext, prior_weight = w, n_curr, n_ext, robust_component_sd, n_robust = 1 / robust_component_sd^2, posterior_weight))
    })
  },
  mc.cores = ncores
))

posterior_weight$title <- "Robust_precision"
posterior_weight$title <- as.factor(posterior_weight$title)

# New facet label names for supp variable
levels(posterior_weight$title) <- c(Robust_precision = latex2exp::TeX("$\\n_{robust}\\, (precision \\, increasing \\, from \\, left \\, to \\, right)$"))

ggplot(data = posterior_weight, aes(x = mu_ext, y = posterior_weight)) +
  geom_line(linewidth = 1.3) +
  facet_nested(. ~ title * n_robust, labeller = labeller(.cols = label_parsed, .rows = label_parsed)) +
  xlab(latex2exp::TeX("$\\bar{y}_{ext}-\\theta$"))+
  ylab("Posterior weight on \n informative component") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 28, family = "serif"),
    axis.title.x = element_text(size = 28, family = "serif"),
    axis.title.y = element_text(size = 30, family = "serif"),
    legend.text = element_text(size = 30, family = "serif"),
    legend.title = element_text(size = 30, family = "serif"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.key.width = unit(11, "line"),
    legend.key.height = unit(1, "line"),
    legend.spacing.x = unit(0, "cm"),
    legend.spacing.y = unit(-0.3, "cm"),
    strip.text.x = element_text(size = 28, family = "serif")
  ) +
  theme(
    panel.border = element_rect(fill = NA, colour = "black"),
    strip.background = element_rect(color = "black", linewidth = 1)
  )
