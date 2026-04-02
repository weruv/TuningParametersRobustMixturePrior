second_deriv <- function(n02, sigma, w_tilda, theta, n, sigma02, y_bar, mu_vague, sigma01, n01, mu_ext) {
  second_deriv <- (1.0 * (n02 / sigma02^2 + n / sigma^2)^(5 / 2) * (1 - w_tilda) * (theta - (n * sigma02^2 * y_bar + mu_vague * n02 * sigma^2) / (n * sigma02^2 + n02 * sigma^2))^2 *
    exp(-(0.5 * (n02 / sigma02^2 + n / sigma^2) * (theta - (n * sigma02^2 * y_bar + mu_vague * n02 * sigma^2) / (n * sigma02^2 + n02 * sigma^2))^2))) / (sqrt(2) * sqrt(pi)) - (1.0 * (n02 / sigma02^2 + n / sigma^2)^(3 / 2) * (1 - w_tilda) *
    exp(-(0.5 * (n02 / sigma02^2 + n / sigma^2) * (theta - (n * sigma02^2 * y_bar + mu_vague * n02 * sigma^2) / (n * sigma02^2 + n02 * sigma^2))^2))) / (sqrt(2) * sqrt(pi)) + (1.0 * (n01 / sigma01^2 + n / sigma^2)^(5 / 2) * (w_tilda) * (theta - (n * sigma01^2 * y_bar + mu_ext * n01 * sigma^2) / (n * sigma01^2 + n01 * sigma^2))^2 *
    exp(-(0.5 * (n01 / sigma01^2 + n / sigma^2) * (theta - (n * sigma01^2 * y_bar + mu_ext * n01 * sigma^2) / (n * sigma01^2 + n01 * sigma^2))^2))) / (sqrt(2) * sqrt(pi)) - (1.0 * (n01 / sigma01^2 + n / sigma^2)^(3 / 2) * (w_tilda) *
    exp(-(0.5 * (n01 / sigma01^2 + n / sigma^2) * (theta - (n * sigma01^2 * y_bar + mu_ext * n01 * sigma^2) / (n * sigma01^2 + n01 * sigma^2))^2))) / (sqrt(2) * sqrt(pi))
  return(data.frame(theta = theta, second_deriv = second_deriv))
}