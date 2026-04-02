TIErate_hybridcontrol <- function(Delta, n_c, n_t, robust_component_location, mu_ext, n_ext, CONTROL.SIGMA, TREAT.SIGMA, ext.SIGMA, n_robust, w, H0Diff, muE, cE, nsim) {
  set.seed(123)
  dat_control <- rnorm(nsim, mean = Delta, sd = CONTROL.SIGMA / sqrt(n_c))
  dat_treat <- rnorm(nsim, mean = Delta + H0Diff, sd = TREAT.SIGMA / sqrt(n_t))

  TIE_rate <- mean(sapply(1:length(dat_control), function(z) {
    if (robust_component_location == "current_observed_mean") {
      prior_control <- mixcombine(mixnorm(informative = c(1, c(mu_ext, ext.SIGMA / sqrt(n_ext))), sigma = CONTROL.SIGMA), mixnorm(vague = c(1, c(dat_control[z], n_robust)), sigma = CONTROL.SIGMA, param = "mn"), weight = c(w, 1 - w))
    } else if (robust_component_location == "ext") {
      prior_control <- mixcombine(mixnorm(informative = c(1, c(mu_ext, ext.SIGMA / sqrt(n_ext))), sigma = CONTROL.SIGMA), mixnorm(vague = c(1, c(mu_ext, n_robust)), sigma = CONTROL.SIGMA, param = "mn"), weight = c(w, 1 - w))
    }
    post.control <- suppressMessages(postmix(prior_control, m = dat_control[z], n = n_c))

    post.treat <- mixnorm(
      comp1 = c(1, dat_treat[z], TREAT.SIGMA / sqrt(n_t)),
      sigma = TREAT.SIGMA
    ) # assumes flat prior
    pmixdiff(
      mix1 = post.treat,
      mix2 = post.control,
      q = muE,
      lower.tail = TRUE
    ) <= 0.025
  }))
   return(TIE_rate)
}
