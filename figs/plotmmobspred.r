# Calculate observed vs predicted, and r2glmm, for the mixed models

library(r2glmm)
library(sjstats)
library(lme4)
library(glmmTMB)
library(purrr)

r2beta(fit_bbs_bcr[[1]]$model)
r2(fit_bbs_bcr[[2]]$model)


# RMSEs of all models


get_rmse <- function(x,nor) map_dbl(x, function(y) rmse(y$model, normalized = nor))
get_rmse(fit_bbs_bcr, TRUE)
get_rmse(fit_bbs_huc, TRUE)
get_rmse(fit_bbs_tnc, TRUE)

# obs vs pred
pred <- predict(fit_bbs_bcr[[2]]$model)
obs <- fit_bbs_bcr[[1]]$model@frame[,1]

# Function to (1) get normalized rmse (2) get observed and predicted values (3) ggplot them with rmse depicted on plot

obs_pred_plot <- function(x, name) {
  require(ggplot2)
  x <- x$model
  rmse_raw <- rmse(x, normalized = FALSE)
  rmse_norm <- rmse(x, normalized = TRUE)
  rmse_label <- paste0(' RMSE: ', round(rmse_norm,2), ' (', round(rmse_raw,2), ' raw)')
  pred <- predict(x)
  if (inherits(x, 'lmerMod')) {
    obs <- x@frame[,1]
  } else {
    obs <- x$frame[,1]
  }
  dat <- data.frame(pred,obs)
  axisrange <- range(c(obs,pred), na.rm = TRUE)
  ggplot(dat, aes(x=pred, y=obs)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = 'indianred') +
    geom_text(data = data.frame(pred = -Inf, obs = Inf, lab = rmse_label), aes(label = lab), hjust = 0, vjust = 1) +
    theme_bw() +
    scale_x_continuous(limits = axisrange, name = 'Predicted') +
    scale_y_continuous(limits = axisrange, name = 'Observed') +
    ggtitle(name)
}

bio_titles <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'beta PD', 'gamma PD', 'alpha FD', 'beta FD', 'gamma FD')

obsplot_bbs_bcr <- map2(fit_bbs_bcr, bio_titles, obs_pred_plot)
obsplot_bbs_huc <- map2(fit_bbs_huc, bio_titles, obs_pred_plot)
obsplot_bbs_tnc <- map2(fit_bbs_tnc, bio_titles, obs_pred_plot)

bio_titles_fia <- c('alpha TD', 'alpha TD abundance', 'alpha PD', 'alpha PD abundance', 'alpha FD', 'alpha FD abundance', 'beta TD', 'beta TD abundance', 'gamma TD', 'gamma TD abundance', 'gamma PD', 'gamma PD abundance', 'gamma FD', 'gamma FD abundance')

obsplot_fia_bcr <- map2(fit_fia_bcr, bio_titles_fia, obs_pred_plot)
obsplot_fia_huc <- map2(fit_fia_huc, bio_titles_fia, obs_pred_plot)
obsplot_fia_tnc <- map2(fit_fia_tnc, bio_titles_fia, obs_pred_plot)

obsplot_fia_huc_incidence <- obsplot_fia_huc[c(1,7,9,3,11,5,13)]
obsplot_fia_huc_abundance <- obsplot_fia_huc[c(2,8,10,4,12,6,14)]
obsplot_fia_bcr_incidence <- obsplot_fia_bcr[c(1,7,9,3,11,5,13)]
obsplot_fia_bcr_abundance <- obsplot_fia_bcr[c(2,8,10,4,12,6,14)]
obsplot_fia_tnc_incidence <- obsplot_fia_tnc[c(1,7,9,3,11,5,13)]
obsplot_fia_tnc_abundance <- obsplot_fia_tnc[c(2,8,10,4,12,6,14)]
lay_mat <- rbind(c(1,2,3),
                 c(4,NA,5),
                 c(6,NA,7))

# Draw figures

library(gridExtra)

fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/observed_predicted_plots/'

png(file.path(fpfig, 'BBS_HUC4_obsvspred.png'), height = 10, width = 10, res = 400, units = 'in')
grid.arrange(grobs = obsplot_bbs_huc, nrow = 3)
dev.off()
png(file.path(fpfig, 'BBS_BCR_obsvspred.png'), height = 10, width = 10, res = 400, units = 'in')
grid.arrange(grobs = obsplot_bbs_bcr, nrow = 3)
dev.off()
png(file.path(fpfig, 'BBS_TNC_obsvspred.png'), height = 10, width = 10, res = 400, units = 'in')
grid.arrange(grobs = obsplot_bbs_tnc, nrow = 3)
dev.off()

png(file.path(fpfig, 'FIA_HUC4_obsvspred_incidence.png'), height = 10, width = 10, res = 400, units = 'in')
grid.arrange(grobs = obsplot_fia_huc_incidence, layout_matrix = lay_mat)
dev.off()
png(file.path(fpfig, 'FIA_BCR_obsvspred_incidence.png'), height = 10, width = 10, res = 400, units = 'in')
grid.arrange(grobs = obsplot_fia_bcr_incidence, layout_matrix = lay_mat)
dev.off()
png(file.path(fpfig, 'FIA_TNC_obsvspred_incidence.png'), height = 10, width = 10, res = 400, units = 'in')
grid.arrange(grobs = obsplot_fia_tnc_incidence, layout_matrix = lay_mat)
dev.off()

png(file.path(fpfig, 'FIA_HUC4_obsvspred_abundance.png'), height = 10, width = 10, res = 400, units = 'in')
grid.arrange(grobs = obsplot_fia_huc_abundance, layout_matrix = lay_mat)
dev.off()
png(file.path(fpfig, 'FIA_BCR_obsvspred_abundance.png'), height = 10, width = 10, res = 400, units = 'in')
grid.arrange(grobs = obsplot_fia_bcr_abundance, layout_matrix = lay_mat)
dev.off()
png(file.path(fpfig, 'FIA_TNC_obsvspred_abundance.png'), height = 10, width = 10, res = 400, units = 'in')
grid.arrange(grobs = obsplot_fia_tnc_abundance, layout_matrix = lay_mat)
dev.off()


# Confidence intervals of coefficients ------------------------------------

load('C:/Users/Q/Dropbox/projects/nasabiodiv/mmcis_bbs.RData')
ci_bbs <- rbind(do.call(rbind, ci_bbs_bcr),
                do.call(rbind, ci_bbs_huc),
                do.call(rbind, ci_bbs_tnc))
ci_bbs$region <- rep(c('BCR','HUC4','TNC'), each = 72)
ci_bbs$resp <- rep(factor(bio_titles, levels = bio_titles), each = 8) # Keep order the same as the maps.

pred_labels <-  c("bio12_5k_100_mean" = 'Precip mean',
                  "bio12_5k_100_sd" = 'Precip stdev',
                  "bio1_5k_100_mean" = 'Temp mean', 
                  "dhi_gpp_5k_100_sd" = 'GPP stdev', 
                  "elevation_5k_100_sd" = 'Elevation stdev', 
                  "geological_age_5k_100_diversity" = 'Geol diversity', 
                  "human_footprint_5k_100_mean" = 'Footprint mean',
                  "soil_type_5k_100_diversity" = 'Soil type diversity')
ci_bbs$pred <- pred_labels[as.character(ci_bbs$pred)]

map(c('BCR','HUC4','TNC'), function(reg) {
p <- ggplot(subset(ci_bbs, region == reg), aes(x = pred, ymin = q025, ymax = q975)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'blue') +
  geom_errorbar(width = 0.2, size = 1.1) +
  facet_wrap(~ resp, scales = 'free_y') +
  theme_bw() +
  scale_y_continuous(name = 'Overall slope') +
  scale_x_discrete(name = 'Predictor variable') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave(file.path(fpfig, paste0('BBS_',reg,'coefficientCIs.png')), p, height = 9, width = 10, dpi = 400)
})


# Confidence intervals of coeffs (FIA) ------------------------------------

bio_titles_fia <- c('alpha TD', 'alpha TD abundance', 'alpha PD', 'alpha PD abundance', 'alpha FD', 'alpha FD abundance', 'beta TD', 'beta TD abundance', 'gamma TD', 'gamma TD abundance', 'gamma PD', 'gamma PD abundance', 'gamma FD', 'gamma FD abundance')

load('C:/Users/Q/Dropbox/projects/nasabiodiv/mmcis_fia.RData')
ci_fia <- rbind(do.call(rbind, ci_fia_bcr),
                do.call(rbind, ci_fia_huc),
                do.call(rbind, ci_fia_tnc))
ci_fia$region <- rep(c('BCR','HUC4','TNC'), each = nrow(ci_fia)/3)
ci_fia$resp <- rep(factor(bio_titles_fia, levels = bio_titles_fia), each = 8) # Keep order the same as the maps.

pred_labels <-  c("bio12_5k_100_mean" = 'Precip mean',
                  "bio12_5k_100_sd" = 'Precip stdev',
                  "bio1_5k_100_mean" = 'Temp mean', 
                  "dhi_gpp_5k_100_sd" = 'GPP stdev', 
                  "elevation_5k_100_sd" = 'Elevation stdev', 
                  "geological_age_5k_100_diversity" = 'Geol diversity', 
                  "human_footprint_5k_100_mean" = 'Footprint mean',
                  "soil_type_5k_100_diversity" = 'Soil type diversity')
ci_fia$pred <- pred_labels[as.character(ci_fia$pred)]
