# Main code below by S. Calderazzo
require(RBesT)
require(parallel)
require(cubature)
require(mvtnorm)
require(VGAM)
library(truncnorm)


tIe <- function(alpha.ev, At, Ac, mupc, mupt, sigma, nt, nc, eval, th0 = 0) {
  1 - pnorm((th0 + mupc * (Ac / (Ac + 1)) - mupt * (At / (At + 1)) - eval * ((1 / (At + 1)) - (1 / (Ac + 1))) +
    qnorm(1 - alpha.ev) * sqrt(((sigma^2 / nc) / (1 + Ac)) + ((sigma^2 / nt) / (1 + At)))) /
    sqrt(((sigma^2 / nc) / (1 + Ac)^2) + ((sigma^2 / nt) / (1 + At)^2)))
}


p.h1.precomp <- function(pow, nc, nt, n0c, n0t, a0t = 0.5, a0c = 0.5, b0t = 0.5, b0c = 0.5, th0 = 0) {
  sapply(0:(pow * n0c + nc), function(xS) {
    sapply(0:(pow * n0t + nt), function(xT) {
      post.treat <- mixbeta(c(1, a0t + xT, b0t + pow * n0t + nt - xT))
      post.control <- mixbeta(c(1, a0c + xS, b0c + pow * n0c + nc - xS))
      out <- try(pmixdiff(post.treat, post.control, th0, lower.tail = FALSE), silent = TRUE)
      if (inherits(out, "try-error")) {
        out <- 1 - fisher.test(matrix(c(xS, pow * n0c + nc - xS, xT, pow * n0t + nt - xT), 2, 2), alternative = "less")$p.value
      }
      as.numeric(out)
    })
  })
}


irtest <- function(dec, oc, c.alpha = 0.975, c.beta = 0.025) {
  c.alpha * (1 - oc) * dec + c.beta * (oc) * (1 - dec)
}

cThreshold <- function(yc, yt, nc, nt, priorPars, outcome = "normal",
                       mc.results = NULL, g = 0.05, sigma = 1, alpha.b = 0.025,
                       alpha.up = 1, alpha.low = 0, a0t = 0.5, a0c = 0.5, b0t = 0.5, b0c = 0.5, th0 = 0,
                       yc.grid = c(seq(-10, -2.5, by = 0.5), seq(-2, 2, by = 0.01), seq(2.5, 10, by = 0.5)),
                       pH1array = NULL, delta1 = NULL, delta2 = NULL, dE) {
  if (outcome == "normal") {
    mupc <- priorPars[1, 1] + dE
    mupt <- priorPars[1, 2]
    sigmapc <- priorPars[2, 1]
    sigmapt <- priorPars[2, 2]

    Ac <- sigma^2 / (nc * sigmapc^2)
    At <- sigma^2 / (nt * sigmapt^2)
    eval <- yc

    # Est. TIE under informative
    tIe.info <- tIe(alpha.ev = alpha.b, At = At, Ac = Ac, mupc = mupc, mupt = mupt, sigma = sigma, nt = nt, nc = nc, eval = eval, th0 = th0)

    kappa.cal <- 1 - pnorm((th0 + mupc * (Ac / (Ac + 1)) - mupt * (At / (At + 1)) - eval * ((1 / (At + 1)) - (1 / (Ac + 1))) +
      qnorm(1 - alpha.b) * sqrt(((sigma^2 / nc) / (1 + Ac)) + ((sigma^2 / nt) / (1 + At)))) * (At + 1) /
      sqrt((sigma^2 / nc) + (sigma^2 / nt)))



    eval1 <- mupc + c(-delta1, delta1)
    eval2 <- mupc + c(-delta2, delta2)

    kappa <- ifelse(eval >= eval1[1] & eval <= eval1[2], kappa.cal,
      ifelse(eval >= eval2[1] & eval <= eval2[2], pmax(pmin(kappa.cal, alpha.up), alpha.low),
        alpha.b
      )
    )
  }


  kappa
}

# adapted from studyPrior package
power_par <- function(dat, n, sigma = 1, priorPars, outcome, bound = FALSE, b0 = 0.5, a0 = 0.5) {
  if (outcome == "normal") {
    sD <- sigma^2 / n
    d <- priorPars[2]^2 / (pmax((dat - priorPars[1])^2, sD + priorPars[2]^2) - sD)
    if (bound) d <- pmin(d, priorPars[2]^2 / (sigma^2 / n))
  }
  if (outcome == "binomial") {
    optfn <- function(datC) {
      lik.d <- function(d) VGAM::dbetabinom.ab(x = datC, size = n, d * (priorPars[1]) + a0, d * (priorPars[2]) + b0) # uniform historical prior

      opd <- BB::spg(
        par = .005,
        fn = lik.d,
        lower = 0,
        upper = 1,
        control = list(
          maximize = TRUE,
          trace = FALSE
        )
      )

      if (opd$convergence != 0) print(opd)

      ds <- opd$par

      return(ds)
    }
    d <- sapply(dat, optfn)
    if (bound) d <- pmin(d, n / (priorPars[1] + priorPars[2]))
  }
  return(d)
}


pH1 <- function(yc, yt, nc, nt, priorPars, w.mix, sigma = 1, prior = "info", bound = FALSE, pH1array = NULL,
                a0t = 0.5, a0c = 0.5, b0t = 0.5, b0c = 0.5, n0 = 1, th0 = 0, d1 = NULL, delta.v = NULL, w = 0.5, robust_comp_location, dE,robust_comp_sd) {
  mupc <- priorPars[1, 1] + dE
  mupt <- priorPars[1, 2]
  sigmapc <- priorPars[2, 1]
  sigmapt <- priorPars[2, 2]

  if (prior == "mix" | prior == "mix2") {
    if (prior == "mix" | prior == "mix2") {
      if (sigmapt < 10^3) {
        mix.treat.unit <- mixnorm(info = c((w.mix), mupt, sigmapt), unit = c(1 - w.mix, mupt, sigma), sigma = sigma)
      } ###
      else {
        mix.treat.unit <- mixnorm(info = c(1, mupt, 10^3), sigma = sigma)
      }

      if (robust_comp_location == "ext") {
        if (sigmapc < 10^3) {
          mix.control.unit <- mixnorm(info = c((w.mix), mupc, sigmapc), unit = c(1 - w.mix, mupc, robust_comp_sd), sigma = sigma)
        } ###
        else {
          mix.control.unit <- mixnorm(info = c(1, mupc, 10^3), sigma = sigma)
        } #

        p.h1 <- sapply(1:length(yc), function(z) {
          post.control.mix <- suppressMessages(postmix(mix.control.unit, m = yc[z], n = nc))
          post.treat <- suppressMessages(postmix(mix.treat.unit, m = yt[z], n = nt))
          pmixdiff(post.treat, post.control.mix, th0, lower.tail = FALSE)
        })
      } else if (robust_comp_location == "current.mean") {
        p.h1 <- sapply(1:length(yc), function(z) {
          if (sigmapc < 10^3) {
            mix.control.unit <- mixnorm(info = c((w.mix), mupc, sigmapc), unit = c(1 - w.mix, yc[z], robust_comp_sd), sigma = sigma)
          } ### 
          else {
            mix.control.unit <- mixnorm(info = c(1, mupc, 10^3), sigma = sigma)
          } # needs more elegant solution
          post.control.mix <- suppressMessages(postmix(mix.control.unit, m = yc[z], n = nc))
          post.treat <- suppressMessages(postmix(mix.treat.unit, m = yt[z], n = nt))
          pmixdiff(post.treat, post.control.mix, th0, lower.tail = FALSE)
        })
      }
    }
    if (prior == "mix2") {
      p.h1 <- sapply(1:length(yc), function(z) {
        if (sigmapc < 10^3) {
          w.mix <- 1 / 2
          Ac <- sigma^2 / (nc * sigma^2)

          D1 <- ((Ac + 1) / Ac) * qnorm(1 - alpha.b) * (sqrt(sigma^2 / nc + sigma^2 / nt) - sqrt((sigma^2 / nc) / (Ac + 1) + sigma^2 / nt))

          Ac1 <- sigma^2 / (nc * sigmapc^2)

          D2 <- ((Ac1 + 1) / Ac1) * (qnorm(1 - alpha.low) * sqrt(sigma^2 / nc + sigma^2 / nt) - qnorm(1 - alpha.b) * sqrt((sigma^2 / nc) / (Ac1 + 1) + sigma^2 / nt))
          D4 <- ((Ac1 + 1) / Ac1) * (qnorm(1 - alpha.up) * sqrt(sigma^2 / nc + sigma^2 / nt) - qnorm(1 - alpha.b) * sqrt((sigma^2 / nc) / (Ac1 + 1) + sigma^2 / nt))

          mix.control.unit <- mixnorm(info = c(w.mix, mupc, sigmapc), delta1 = c(w.mix, yc[z] + D1, sigma), sigma = sigma)
        } else {
          mix.control.unit <- mixnorm(info = c(1, mupc, 10^3), sigma = sigma)
        } # needs more elegant solution

        if (sigmapt < 10^3) {
          mix.treat.unit <- mixnorm(info = c((1 - w.mix), mupt, sigmapt), unit = c(w.mix, mupt, sigma), sigma = sigma)
        } else {
          mix.treat.unit <- mixnorm(info = c(1, mupt, 10^3), sigma = sigma)
        }

        post.control.mix <- suppressMessages(postmix(mix.control.unit, m = yc[z], n = nc))
        post.treat <- suppressMessages(postmix(mix.treat.unit, m = yt[z], n = nt))
        pmixdiff(post.treat, post.control.mix, th0, lower.tail = FALSE)
      })
    }
  } else {
    if (prior == "EB") {
      if (sigmapc < 10^3) {
        dc <- power_par(dat = yc, n = nc, sigma = sigma, priorPars = c(mupc, sigmapc), outcome = "normal", bound = bound)
        Ac <- dc * sigma^2 / (nc * sigmapc^2)
      } else {
        Ac <- 0
      }

      if (sigmapt < 10^3) {
        dt <- power_par(dat = yt, n = nt, sigma = sigma, priorPars = c(mupt, sigmapt), outcome = "normal", bound = bound)
        At <- dt * sigma^2 / (nc * sigmapt^2)
      } else {
        At <- 0
      }
    } else {
      Ac <- sigma^2 / (nc * sigmapc^2)
      At <- sigma^2 / (nt * sigmapt^2)
    }

    pthr <- ifelse(rep(sigmapt, length(yc)) > 0,
      ((mupt * At + yt) / (At + 1) - (mupc * Ac + yc) / (Ac + 1)) / sqrt(((sigma^2 / nc) / (1 + Ac)) + ((sigma^2 / nt) / (1 + At))),
      (mupt - (mupc * Ac + yc) / (Ac + 1)) / sqrt(((sigma^2 / nc) / (1 + Ac)))
    )

    p.h1 <- pnorm(th0, pthr, 1, lower.tail = FALSE)
  }
  if (prior == "SampUnif") {
    mupt.post <- (mupt * At + yt) / (At + 1)
    mupc.post <- (mupc * Ac + yc) / (Ac + 1)
    sigmapc2.post <- ((sigma^2 / nc) / (1 + Ac))
    sigmapt2.post <- ((sigma^2 / nt) / (1 + At))

    m1 <- mupc.post
    m2 <- mupt.post
    t12 <- sigmapc2.post
    t22 <- sigmapt2.post
    c1 <- (-1)
    c2 <- 1

    th1 <- c1 * m1 + c2 * m2
    s12 <- c1^2 * t12 + c2^2 * t22
    th2 <- m2
    s22 <- t22

    m1 <- (th0 - th1) / sqrt(s12)
    omega <- (mupc - d1 - th2) / sqrt(s22)
    delta <- (mupc + d1 - th2) / sqrt(s22)
    rho <- c2 * sqrt(s22) / sqrt(s12)
    cmat <- matrix(c(1, rho, rho, 1), 2, 2)

    p.h1 <- sapply(1:length(yt), function(s) {
      1 - (pmvnorm(lower = -Inf, upper = c(m1[s], delta[s]), corr = cmat) - pmvnorm(lower = c(-Inf, -Inf), upper = c(m1[s], omega[s]), corr = cmat)) /
        (pnorm(delta[s]) - pnorm(omega[s]))
    })
  }
  if (prior == "SampPM") {
    weights <- c(w, 1 - w)

    Sigma <- matrix(c(sigmapc^2 + (sigma^2 / nc), -(sigma^2 / nc), -(sigma^2 / nc), sigma^2 / nc + sigma^2 / nt), 2, 2)
    # comp1
    marg.1 <- dmvnorm(cbind(yc, yt - yc), mean = c(mupc, delta.v[1]), sigma = Sigma, log = TRUE, checkSymmetry = FALSE)
    # comp2
    marg.2 <- dmvnorm(cbind(yc, yt - yc), mean = c(mupc, delta.v[2]), sigma = Sigma, log = TRUE, checkSymmetry = FALSE)

    margT <- exp(log(weights[1]) + marg.1) + exp(log(weights[2]) + marg.2)
    marg1 <- exp(log(weights[1]) + marg.1)
    post.weights <- cbind(exp(log(marg1) - log(margT)), 1 - exp(log(marg1) - log(margT)))

    p.h1 <- post.weights[, 2]
  }
  p.h1
}
