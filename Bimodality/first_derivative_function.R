first_derivative <- function(n, sigma, n0_1, sigma0_1, theta, m1, y_bar) {
  derivative <- (-1) * (((1 / sqrt(2 * pi)) * ((sqrt(n / sigma^2 + n0_1 / sigma0_1^2)) * (n / sigma^2 + n0_1 / sigma0_1^2) * (theta - ((n0_1 * m1 * sigma^2 + n * y_bar * sigma0_1^2) / (n0_1 * sigma^2 + n * sigma0_1^2))) * exp(((-0.5) * (n / sigma^2 + n0_1 / sigma0_1^2)) * (theta - ((n0_1 * m1 * sigma^2 + n * y_bar * sigma0_1^2) / (n0_1 * sigma^2 + n * sigma0_1^2)))^2))
  ))
  return(data.frame(theta, derivative))
}