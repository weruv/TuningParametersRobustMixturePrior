posterior_dens_func <- function(n0, n1, sigma, sigma0, y_bar, m1, theta) {
  dens <- (1 / sqrt(2 * pi)) * ((sqrt(n1 / sigma^2 + n0 / sigma0^2)) * exp(((-0.5) * (n1 / sigma^2 + n0 / sigma0^2)) * (theta - ((n0 * m1 * sigma^2 + n1 * y_bar * sigma0^2) / (n0 * sigma^2 + n1 * sigma0^2)))^2))
  return(dens)
}