library(RBesT)
library(parallel)
library(latex2exp)
library(tidyr)

# ==============================================================================================================
## OBM computation as in the paper
# ==============================================================================================================

source("posterior_dens_func.R")                                 # posterior density
source("first_derivative_function.R")                           # computes first derivative
source("second_derivative_function.R")                          # computes second derivative
source("importance_of_modes_bimodality.R")                      # computes OBM

mu_ext <- seq(0, 2, by = 0.01)                                 ## prior-data conflict, x-axis in Figure 2
current_mu <- 0                                                 # current data true mean
robust_component_sd <- c(sqrt(0.5), 1, 2, 5, 20)                ## dispersion of robust component
n_ext <- c(15)                                                 # external data sample size
n_curr <- c(20)                                                 # current sample size
w <- seq(0, 1, by = 0.01)                                       # weight given to the informative component (0 < weight < 1).
SIGMA <- 1                                                      # current sigma, assumed as known...
sigma <- 1                                                      # historical sigma, assumed as known...

robust_comp_location <- c("ext", "null hypothesis")            ## location of robust component

grid <- expand.grid(w = w, mu_ext = mu_ext, robust_component_sd = robust_component_sd, robust_comp_location = robust_comp_location, n_ext = n_ext, n_curr = n_curr)

ncores <- 64

OBM_results <- do.call(rbind, mclapply( split(grid, seq_len(nrow(grid))),
                                        function(params)  {
        with(params, {                                     
        RhpcBLASctl::blas_set_num_threads(1)
          
        if (robust_comp_location == "ext") {
        mu_robust <- mu_ext
        } else if (robust_comp_location == "null hypothesis") {
         mu_robust <- 0 # null hypothesis as specified in the paper
        }
        results_OBM <- importance_of_modes(w = w, ybar = current_mu, sigma_info = sigma / sqrt(n_ext), sigma_robust = robust_component_sd, mu_ext = mu_ext, mu_robust = mu_robust, sigma = SIGMA, ext_sd = sigma, n = n_curr)
        results_OBM <- cbind(results_OBM, robust_comp_location = robust_comp_location, n_ext = n_ext, n_curr = n_curr)
        return(results_OBM)
        })
      }, mc.cores = ncores))

## ---------------------------------------------------------------------------------------------------------
## robust_component_sd in terms of number of samples
OBM_results$n_robust <- 1 / OBM_results$sigma_robust^2
## ---------------------------------------------------------------------------------------------------------

## ---------------------------------------------------------------------------------------------------------
# plot (Figure 2)
## subset to plot for unit-info robust component
## without the subset, you obtain Figure A.9
## ---------------------------------------------------------------------------------------------------------
OBM_results <- subset(OBM_results, n_robust == 1)

## group OBM into several classes
OBM_results$discrete_OBM <- cut(OBM_results$OBM, breaks = seq(from = min(OBM_results$OBM, na.rm = T), to = max(OBM_results$OBM, na.rm = T), length.out = 8), include.lowest = T)
levels(OBM_results$discrete_OBM) <- c(levels(OBM_results$discrete_OBM), "Unimodal")
OBM_results$discrete_OBM[is.na(OBM_results$discrete_OBM)] <- "Unimodal"
OBM_results$discrete_OBM <- relevel(OBM_results$discrete_OBM, ref = "Unimodal")

##
titles <- c(TeX("$\\mu_{robust}=\\bar{y}_{ext}$"), TeX("$\\mu_{robust}=\\bar{y}$"))
conditions <- unique(OBM_results$robust_component_location)

Fig_bimodality_intensity <- lapply(1:length(conditions), function(indexs) {
  p <- OBM_results %>%
    dplyr::filter(robust_component_location == conditions[indexs]) %>%
    ggplot(., aes(x = mu_ext, y = w, fill = discrete_OBM)) +
    geom_raster() +
    ylab("Prior Weight (w) \n on external data") +
    xlab(TeX("$\\bar{y}_{ext}-\\theta$")) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 28, family = "serif"),
      axis.title.x = element_text(size = 28, family = "serif"),
      axis.title.y = element_text(size = 30, family = "serif"),
      legend.text = element_text(size = 30, family = "serif"),
      legend.title = element_text(size = 30, family = "serif"),
      legend.direction = "horizontal", legend.key.size = unit(2, "cm"), legend.key.height = unit(1, "lines"),
      legend.position = "bottom",
      strip.text.x = element_text(size = 18),
      strip.text.y = element_text(size = 21),
      plot.title = element_text(hjust = 0.5, size = 28)
    ) +
    theme(
      axis.text = element_text(color = c("black", "black", "black", "black", "black")),
      axis.ticks = element_line(color = c("black", "black", "black", "black", "black"))
    ) +
    labs(fill = "O'Hagan's Bimodality \n metric") +
    ggtitle(titles[indexs]) +
    scale_x_continuous(breaks = c(0, 1, 2)) +
    scale_fill_brewer(palette = "RdYlBu", na.value = "#5b82b4", direction = -1) +
    theme(panel.background = element_rect(fill = "white", colour = "black")) +
    guides(fill = "none")
})
