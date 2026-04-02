# ============================================================
# code below by S. Weber
# ============================================================
# Student-t prior + Normal likelihood (known sigma)
# Normal-mixture approximation as RBesT::mixnorm objects
#   - Methods: "gl_prob" (Gauss–Laguerre, Gamma prob-weight) and "quantile"
#   - Routes:  "direct" (λ reweight by p(m|λ)) and "indirect" (prior -> update + reweight)
# Notation: m0 = prior mean, s0 = prior scale (sd), nu = df
# ============================================================

# Optional dependency for Gauss–Laguerre ("gl_prob"):
# install.packages("statmod", repos = "https://cloud.r-project.org")
# RBesT (for mixnorm):
# install.packages("RBesT", repos = "https://cloud.r-project.org")

suppressPackageStartupMessages({
  have_statmod <- requireNamespace("statmod", quietly = TRUE)
  have_RBesT   <- requireNamespace("RBesT",   quietly = TRUE)
})

if (!have_RBesT) stop("This script requires RBesT. Please install with install.packages('RBesT').")

# -------------------------------
# 1) λ-node generators (Γ(shape=nu/2, rate=nu/2))
# -------------------------------
lambda_nodes_gl_prob <- function(nu, nodes, normalize = TRUE) {
  if (!have_statmod) stop('GL-prob requires the "statmod" package.')
  alpha <- nu/2; beta <- nu/2
  quad  <- statmod::gauss.quad.prob(nodes, dist = "gamma", alpha = alpha, beta = 1)
  lam   <- quad$nodes / beta
  w     <- quad$weights
  if (normalize) w <- w / sum(w)
  list(lambda = lam, weight = w, method = "gl_prob")
}

lambda_nodes_quantile <- function(nu, nodes, normalize = TRUE) {
  alpha <- nu/2; beta <- nu/2
  u   <- (seq_len(nodes) - 0.5) / nodes
  lam <- qgamma(u, shape = alpha, rate = beta)
  w   <- rep(1/nodes, nodes)
  if (normalize) w <- w / sum(w)
  list(lambda = lam, weight = w, method = "quantile")
}

# -------------------------------
# 2) PRIOR components from λ (θ|λ ~ N(m0, s0^2/λ))
# -------------------------------
prior_components_from_lambda <- function(nu, m0, s0, lambda, base_w, normalize = TRUE) {
  stopifnot(nu > 0, s0 > 0, length(lambda) == length(base_w))
  w <- if (normalize) base_w / sum(base_w) else base_w
  data.frame(
    weight = as.numeric(w),
    mean   = rep(m0, length(lambda)),
    sd     = as.numeric(s0 / sqrt(lambda)),
    lambda = as.numeric(lambda),
    route  = "prior"
  )
}

# -------------------------------
# 3) Convert components -> RBesT::mixnorm
#    mixnorm expects triplets (w, m, s) and a reference scale `sigma`; we use ms-param.
# -------------------------------
components_to_mixnorm <- function(comps, sigma_ref) {
  stopifnot(requireNamespace("RBesT", quietly = TRUE))
  stopifnot(all(c("weight","mean","sd") %in% names(comps)))
  K <- nrow(comps)
  # build list of named triplets
  comp_list <- setNames(lapply(seq_len(K), function(i) {
    c(comps$weight[i], comps$mean[i], comps$sd[i])
  }), paste0("c", seq_len(K)))
  # call mixnorm(..., sigma = sigma_ref, param = "ms")
  do.call(RBesT::mixnorm, c(comp_list, list(sigma = sigma_ref, param = "ms")))
  # mixnorm docs & interface details: weights & (m,s) triplets; requires sigma reference scale. [1](https://opensource.nibr.com/RBesT/reference/mixnorm.html)[2](https://rdrr.io/cran/RBesT/man/mixnorm.html)
}

# -------------------------------
# 4) Top-level builders returning RBesT::mixnorm
# -------------------------------
# Prior (n=0)
mixnorm_t_prior <- function(nu, m0, s0, nodes = 40,
                            method = c("gl_prob","quantile"),
                            sigma_ref = 1, normalize = TRUE) {
  method <- match.arg(method)
  lw <- if (method == "gl_prob") lambda_nodes_gl_prob(nu, nodes) else lambda_nodes_quantile(nu, nodes)
  comps <- prior_components_from_lambda(nu, m0, s0, lw$lambda, lw$weight, normalize = normalize)
  components_to_mixnorm(comps, sigma_ref = sigma_ref)
}
