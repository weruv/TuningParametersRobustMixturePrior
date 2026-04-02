library(RBesT)
library(parallel)
library(ggplot2)
library(latex2exp)
library(ggh4x)
library(cowplot)

## ===========================================================================================================
## average TIE rate and power- hybrid control
## ===========================================================================================================
source("function_averageOCs.R")

nt <- 20                                             # sample size in treatment arm
nc <- 20                                             # sample size in control arm
n_ext <- 15                                         # historical data sample size
mu_ext <- seq(-3, 3, by = 0.01)                     ## prior data conflict in paper, x-axis of Figure 1

design_prior_type <- c("info", "mix", "unit")         ## different kinds of design priors

grid <- expand.grid(nt = nt, nc = nc, n_ext = n_ext, mu_ext = mu_ext, design_prior_type = design_prior_type)

ncores <- 64

average_OC <- do.call(rbind, mclapply(split(grid, seq_len(nrow(grid))),
  function(params) {
    with(params, {
      sigma <- 1
      sig2 <- sigma^2                       # variance of one observation

      sigmapc <- sqrt(sig2 / n_ext)         # SD control arm prior
      sigmapt <- sqrt(sig2 / 0)             # SD treat arm prior
      muc <- 0                              # prior mean control
      mut <- 0.4                            # prior mean treatment
      priorPars <- rbind(c(muc, mut), c(sigmapc, sigmapt))


      alpha.b <- 0.025                      # frequentist TIE
      th0 <- 0

      bound <- FALSE #
      nmc <- 1e6                            # number of Monte Carlo samples

      th1 <- 0.8330874                      # point for Power evaluation
      delta <- c(th0, th1)

      grid_mixture <- expand.grid(robust_component_sd = c(1), robust_comp_location = c("ext", "current.mean"), prior_weight = c(0.25, 0.5, 0.75))
      
      set.seed(123)
      Avg.OC <- sapply(delta, function(x) {
        if (design_prior_type == "info") {
          thc <- rnorm(nmc, priorPars[1, 1], priorPars[2, 1])
        } else if (design_prior_type == "mix") {
          thc <- rmix(mixnorm(informative = c(0.5, c(priorPars[1, 1], priorPars[2, 1])), vague = c(1 - 0.5, c(priorPars[1, 1], robust_component_sd)), sigma = 1), nmc)
        } else if (design_prior_type == "unit") {
          thc <- rnorm(nmc, priorPars[1, 1], sqrt(sig2))
        }

        # sample data
        yc <- rnorm(nmc, thc, sigma / sqrt(nc))
        yt <- rnorm(nmc, thc + x, sigma / sqrt(nt))

        # Weight vector, 1 for monte carlo
        wmatT <- rep(1 / nmc, nmc)

        priorPars.v <- priorPars
        priorPars.v[2, ] <- c(Inf, Inf)
        p.h1.unif <- pH1(
          yc = yc, yt = yt, nc = nc, nt = nt,
          priorPars = priorPars.v, w.mix = w.mix, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext
        )

        p.h1.info <- pH1(
          yc = yc, yt = yt, nc = nc, nt = nt,
          priorPars = priorPars, w.mix = w.mix, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext
        )

        p.h1.mix <- mapply(
          function(prior_weight, robust_comp_location, robust_component_sd) {
            pH1(
              prior = "mix", yc = yc, yt = yt, nc = nc, nt = nt,
              priorPars = priorPars, w.mix = prior_weight, th0 = th0, robust_comp_location = robust_comp_location, dE = mu_ext, robust_comp_sd = robust_component_sd
            )
          },
          grid_mixture$prior_weight, grid_mixture$robust_comp_location, grid_mixture$robust_component_sd
        )


        # test decision
        p.h1.list <- cbind(p.h1.info, p.h1.unif, p.h1.mix)

        dec.list <- cbind(
          p.h1.info > (1 - alpha.b),
          p.h1.unif > (1 - alpha.b),
          p.h1.mix > (1 - alpha.b)
        )


        AvgOC <- wmatT %*% dec.list
        AvgOC
      })
      names_mix <- do.call(paste, c(list("mix"), grid_mixture, sep = "_"))
      rownames(Avg.OC) <- c("info", "freq", names_mix)
      colnames(Avg.OC) <- c("average_TIE", "average_Power")
      
      Avg.OC <- as.data.frame(Avg.OC)
      Avg.OC <- cbind(Avg.OC, mu_ext = mu_ext, nt, nc, n_ext, design_prior_type = design_prior_type)
      Avg.OC$analysis_prior_type <- rownames(Avg.OC)
      
      Avg.OC$robust_location <- ifelse(grepl("mix", Avg.OC$analysis_prior_type), sub(".*_([^_]+)_[^_]+$", "\\1", Avg.OC$analysis_prior_type), NA)
      Avg.OC$robust_sd <- as.numeric(ifelse(grepl("mix", Avg.OC$analysis_prior_type), sub("^[^_]+_([^_]+)_.*", "\\1", Avg.OC$analysis_prior_type), NA))
      Avg.OC$w.mix <- as.numeric(ifelse(grepl("mix", Avg.OC$analysis_prior_type), sub(".*_", "", Avg.OC$analysis_prior_type), NA))
      return(Avg.OC)
    })
  },
  mc.cores = ncores
))

#-------------------------------------------------------------------------------------
## plot (Figure 8 in the paper)
## this plots average TIE rate. To plot average power, change y-axis variable to average_Power
#-------------------------------------------------------------------------------------
average_OC$design_prior_analysis_prior <- paste(average_OC$design_prior_type, average_OC$analysis_prior_type, sep = "_")

# subset the data to plot for weight =0.5
subset_df <- average_OC[!(average_OC$w.mix %in% c(0.25, 0.75) & average_OC$analysis_prior_type == "mix"), ]

Avg_tierate_Fig8a <- ggplot() +
  geom_hline(yintercept = c(0.025), size = 1.2) +
  geom_line(data = subset_df[subset_df$robust_comp_location == "ext" & (subset_df$design_prior_analysis_prior %in% c("info_info", "info_mix", "mix_mix", "unit_mix", "mix_info")), ], aes(x = 0 - mu_ext, y = average_TIE, color = as.factor(design_prior_analysis_prior)), size = 1.2) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 0.08)) +
  scale_color_discrete(labels = c("info_info" = "Informative_Informative", "info_mix" = "Informative_RMP", "mix_mix" = "RMP_RMP", "unit_mix" = "UnitInfo_RMP", "mix_info" = "RMP_Informative")) +
  xlab(TeX("$\\theta_c-\\bar{y}_{ext}$")) +
  ylab("Average TIE rate") +
  theme_bw() +
  ggtitle(TeX("$\\mu_{robust}=\\bar{y}_{ext}$")) +
  theme(
    axis.text = element_text(color = c("black", "black", "black", "black", "black")),
    axis.ticks = element_line(color = c("black", "black", "black", "black", "black"))
  ) +
  theme(
    axis.text = element_text(size = 28, family = "serif"),
    axis.title.x = element_text(size = 28, family = "serif"),
    axis.title.y = element_text(size = 30, family = "serif"),
    legend.text = element_text(size = 30, family = "serif"),
    legend.title = element_text(size = 28, family = "serif"),
    legend.position = c(0.2, 0.85),
    legend.box = "vertical",
    legend.key.width = unit(6, "line"),
    legend.spacing.x = unit(1.2, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    plot.title = element_text(
      hjust = 0.5, # Centers the title (0 = left, 0.5 = center, 1 = right)
      size = 30, # Increases title size
      face = "bold", # (Optional) makes it bold
      margin = margin(b = 20)
    )
  ) +
  guides(color = "none") +
  theme(panel.background = element_rect(fill = "white", colour = NA)) +
  theme(
    axis.line = element_line(colour = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

Avg_tierate_Fig8b <- ggplot() +
  geom_hline(yintercept = c(0.025), size = 1.2) +
  geom_line(data = subset_df[subset_df$robust_comp_location == "current_mean" & (subset_df$design_prior_analysis_prior %in% c("info_info", "info_mix", "mix_mix", "unit_mix", "mix_info")), ], aes(x = 0 - mu_ext, y = average_TIE, color = as.factor(design_prior_analysis_prior)), size = 1.2) +
  scale_color_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 0.08)) +
  scale_color_discrete(labels = c("info_info" = "Informative_Informative", "info_mix" = "Informative_RMP", "mix_mix" = "RMP_RMP", "unit_mix" = "UnitInfo_RMP", "mix_info" = "RMP_Informative")) +
  xlab(TeX("$\\theta_c-\\bar{y}_{ext}$")) +
  ylab("Averate TIE rate") +
  theme_bw() +
  ggtitle(TeX("$\\mu_{robust}=\\bar{y}_{c}$")) +
  theme(
    axis.text = element_text(color = c("black", "black", "black", "black", "black")),
    axis.ticks = element_line(color = c("black", "black", "black", "black", "black"))
  ) +
  theme(
    axis.text = element_text(size = 28, family = "serif"),
    axis.title.x = element_text(size = 28, family = "serif"),
    axis.title.y = element_text(size = 30, family = "serif"),
    legend.text = element_text(size = 30, family = "serif"),
    legend.title = element_text(size = 28, family = "serif"),
    legend.position = c(0.2, 0.85),
    legend.box = "vertical",
    legend.key.width = unit(6, "line"),
    legend.spacing.x = unit(1.2, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
    plot.title = element_text(
      hjust = 0.5, # Centers the title (0 = left, 0.5 = center, 1 = right)
      size = 30, # Increases title size
      face = "bold", # (Optional) makes it bold
      margin = margin(b = 20)
    )
  ) +
  guides(color = "none") +
  theme(panel.background = element_rect(fill = "white", colour = NA)) +
  theme(
    axis.line = element_line(colour = "black", size = 1),
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )

combined_plot <- cowplot::plot_grid(Avg_tierate_Fig8a, NULL, Avg_tierate_Fig8b, nrow = 1, rel_widths = c(1, 0.06, 1), labels = c("(a)", "", "(b)"), align = "hv", axis = "bt", label_size = 20, scale = 0.9, hjust = -1) + panel_border(remove = TRUE)

# extract legend from plot2
legend <- get_legend(
  Avg_tierate_Fig8a +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(title = "Design prior_Analysis prior", nrow = 2))
)
# Combine combined plot and legend using plot_grid()
plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, .1)) + panel_border(remove = TRUE)
