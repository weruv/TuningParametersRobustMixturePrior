importance_of_modes <- function(w, ybar, sigma_info, sigma_robust, mu_ext, mu_robust, sigma, ext_sd, n) {
  ## -----------------------------------------------------
  ## posterior density
  ## -----------------------------------------------------
  priormix <- mixnorm(informative = c(w, c(mu_ext, sigma_info)), vague = c(1 - w, c(mu_robust, sigma_robust)), sigma = ext_sd)
  suppressMessages({
    posteriormix <- postmix(priormix, m = ybar, n = n)
  })
  posterior_dens <- do.call(rbind, lapply(seq(-2, 2, by = 0.01), function(theta) {
    mylist <- list()
    for (i in 1:ncol(posteriormix)) {
      mylist[[i]] <- posteriormix["w", i] * posterior_dens_func(n0 = 1, n1 = n, sigma = sigma, sigma0 = priormix["s", i], y_bar = ybar, m1 = priormix["m", i], theta = theta)
    }
    return(data.frame(theta = theta, post_dens = do.call(sum, mylist)))
  }))

  ## -----------------------------------------------------
  ## first derivative of the posterior
  ## obtain the points where the first derivative changes signs i.e. positive to negative or negative to positive
  ## -----------------------------------------------------
  first_deriv <- do.call(rbind, lapply(seq(-2, 2, by = 0.01), function(theta) {
    mylist <- list()
    for (i in 1:ncol(priormix)) {
      mylist[[i]] <- posteriormix["w", i] * first_derivative(n = n, n0_1 = 1, sigma = 1, sigma0_1 = priormix["s", i], theta = theta, m1 = priormix["m", i], y_bar = ybar)$derivative
    }
    return(data.frame(theta = theta, derivat = do.call(sum, mylist)))
  }))
  ## turning points
  turning_points <- first_deriv$theta[which(diff(sign(first_deriv$derivat)) != 0)]
  results <- data.frame(turning_points = turning_points)

  ## -----------------------------------------------------
  ## second derivative of the posterior
  ## this is to confirm the nature of the points obtained above
  ## if the 2nd derivative < 0, local maximum/mode
  ## if the 2nd derivative > 0, local minimum/antimode
  ## -----------------------------------------------------
  second_derivative <- do.call(rbind, lapply(seq(-2, 2, by = 0.01), function(theta) {
    second_der <- second_deriv(n02 = 1, sigma = sigma, w_tilda = posteriormix["w", 1], theta = theta, n = n, sigma02 = sigma_robust, sigma01 = sigma_info, y_bar = ybar, mu_vague = mu_robust, n01 = 1, mu_ext = mu_ext)
    return(data.frame(theta = theta, second_derivat = second_der$second_deriv))
  }))
  results$second_derivative <- second_derivative$second_derivat[second_derivative$theta %in% results$turning_points]
  results$mode_status <- ifelse(results$second_derivative < 0, "mode", ifelse(results$second_derivative > 0, "antimode", "inflexion"))
  results$posterior <- posterior_dens$post_dens[posterior_dens$theta %in% results$turning_points]

  ## -----------------------------------------------------
  ## compute OBM as in the paper
  ## -----------------------------------------------------
  OBM <- NA
  ratio_of_heights <- NA
  if (nrow(results) > 3) {
    ratio_of_heights <- "multimodal"
    OBM <- "multimodal"
  } else if (nrow(results) == 3) {
    ratio_of_heights <- results$posterior[1] / results$posterior[3]
    OBM <- min(results$posterior[1] / results$posterior[2], results$posterior[3] / results$posterior[2])
  } else {
    ratio_of_heights <- NA
    OBM <- NA
  }
  return(data.frame(mu_ext = mu_ext, mu_robust = mu_robust, sigma_info = sigma_info, sigma_robust = sigma_robust, n = n, w = w, ratio_of_heights = ratio_of_heights, OBM = OBM))
}
