# Code below by S. Calderazzo with minor adjustments by V. Weru
pH1.onearm <- function(yt, nt, priorPars, w.mix, sigma = 1,
                       bound = FALSE, pH1array = NULL, a0t = 0.5, a0c = 0.5, n0 = 1, th0 = 0, sigma_robust, robust_location) {
  mupt <- priorPars[1]
  sigmapt <- priorPars[2]

  if (robust_location == "ext") {
    priorPars.mix <- cbind(priorPars, c(priorPars[1], sigma_robust))
  } else if (robust_location == "H0") {
    priorPars.mix <- cbind(priorPars, c(th0, sigma_robust))
  } else if (robust_location == "current_data_mean") {
    priorPars.mix <- cbind(priorPars, c(yt, sigma_robust))
  }

  weights <- c(w.mix, 1 - w.mix)

  post.var <- 1 / (1 / priorPars.mix[2, ]^2 + nt / sigma^2)
  post.mean <- matrix(NA, length(yt), 2)
  gamma <- (sigma^2 / nt) / ((sigma^2 / nt) + priorPars.mix[2, ]^2)
  post.mean[, 1] <- gamma[1] * priorPars.mix[1, 1] + (1 - gamma[1]) * yt
  post.mean[, 2] <- gamma[2] * priorPars.mix[1, 2] + (1 - gamma[2]) * yt
  dataParPred <- matrix(NA, length(yt), 2)
  dataParPred[, 1] <- sqrt(priorPars.mix[2, 1]^2 + sigma^2 / nt)
  dataParPred[, 2] <- sqrt(priorPars.mix[2, 2]^2 + sigma^2 / nt)

  margT <- exp(log(weights[1]) + dnorm(yt, priorPars.mix[1, 1], dataParPred[, 1], log = TRUE)) + exp(log(weights[2]) + dnorm(yt, priorPars.mix[1, 2], dataParPred[, 2], log = TRUE))
  marg1 <- exp(log(weights[1]) + dnorm(yt, priorPars.mix[1, 1], dataParPred[, 1], log = TRUE))
  post.weights <- cbind(exp(log(marg1) - log(margT)), 1 - exp(log(marg1) - log(margT)))
  wT <- post.weights[, 1]
  p.h1 <- 1 - ((1 - post.weights[, 2]) * pnorm(th0, post.mean[, 1], rep(sqrt(post.var[1]), length(yt))) + post.weights[, 2] * pnorm(th0, post.mean[, 2], rep(sqrt(post.var[2]), length(yt))))
  p.h1
}
