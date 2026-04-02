# analytical power calculation assuming equal sample size in both groups

Powerwo_2sample <- function(treat.mu, control.mu, treat.sigma, control.sigma, n, alpha) {
  power <- pnorm(sqrt(n) * (treat.mu - control.mu) / sqrt(treat.sigma^2 + control.sigma^2) - qnorm(1 - alpha))
}

Powerwo_1sample <- function(current.mu, H0MU, current.sigma, current.n, alpha) {
  power <- pnorm(sqrt(current.n) * (current.mu - H0MU) / current.sigma - qnorm(1 - alpha))
}
