# Fixed effect coefficient plots from spatial mixed models
# QDR/NASABioxgeo/27 Apr 2018

# Edit 02 May: add RMSE and observed vs predicted plots

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv' # Local

coef_bbs <- read.csv(file.path(fp, 'spatial_coef_bbs.csv'), stringsAsFactors = FALSE)
coef_fia <- read.csv(file.path(fp, 'spatial_coef_fia.csv'), stringsAsFactors = FALSE)
spatial_r2s <- read.csv(file.path(fp, 'spatial_r2s.csv'), stringsAsFactors = FALSE)
pred_bbs <- read.csv(file.path(fp, 'spatial_pred_bbs.csv'), stringsAsFactors = FALSE)
pred_fia <- read.csv(file.path(fp, 'spatial_pred_fia.csv'), stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)

# Calculation of RMSE
rmse_bbs <- pred_bbs %>%
  group_by(rv, ecoregion) %>%
  summarize(RMSE = sqrt(mean((observed-Estimate)^2)),
            range_obs = diff(range(observed)),
            rRMSE = RMSE/range_obs)
rmse_fia <- pred_fia %>%
  group_by(rv, ecoregion) %>%
  summarize(RMSE = sqrt(mean((observed-Estimate)^2)),
            range_obs = diff(range(observed)),
            rRMSE = RMSE/range_obs)

rmse_all <- rbind(data.frame(taxon = 'fia', rmse_fia), data.frame(taxon = 'bbs', rmse_bbs))

prednames <- c('elevation_5k_100_sd', 'bio1_5k_100_mean', 'geological_age_5k_100_diversity', 'soil_type_5k_100_diversity', 'bio12_5k_100_mean', 'bio12_5k_100_sd', 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean')
geo_names <- c('elevation_sd','temperature_mean','geol_age_diversity','soil_diversity','precip_mean','precip_sd','gpp_sd','footprint_mean')

bio_titles <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'beta PD', 'gamma PD', 'alpha FD', 'beta FD', 'gamma FD')
bio_names <- c("alpha_richness", "beta_td_sorensen_pa", "gamma_richness",
               "alpha_phy_pa", "beta_phy_pa", "gamma_phy_pa", 
               "alpha_func_pa", "beta_func_pa", "gamma_func_pa")

fia_bio_names_incid <- c("alpha_richness", "beta_td_sorensen_pa", "gamma_richness", 
                         "alpha_phy_pa", "beta_phy_pa", "gamma_phy_pa",
                         "alpha_func_pa", "beta_func_pa", "gamma_func_pa")
fia_bio_names_abund <- c("alpha_effspn", "beta_td_sorensen", "gamma_effspn",
                         "alpha_phy", "beta_phy", "gamma_phy",
                         "alpha_func", "beta_func", "gamma_func")

bio_titles_incidence <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'beta PD', 'gamma PD', 'alpha FD', 'beta FD', 'gamma FD')
bio_titles_abundance <- paste(bio_titles_incidence, 'abundance')

coef_bbs_fixed <- coef_bbs %>% 
  filter(effect == 'fixed', !parameter %in% 'Intercept') %>%
  mutate(predictor = geo_names[match(parameter, prednames)],
         response = factor(bio_titles[match(rv, bio_names)], levels = bio_titles)) 
coef_fia_fixed_incid <- coef_fia %>%
  filter(effect == 'fixed', !parameter %in% 'Intercept', rv %in% fia_bio_names_incid) %>%
  mutate(predictor = geo_names[match(parameter, prednames)],
         response = factor(bio_titles_incidence[match(rv, fia_bio_names_incid)], levels = bio_titles_incidence)) 
coef_fia_fixed_abund <- coef_fia %>%
  filter(effect == 'fixed', !parameter %in% 'Intercept', rv %in% fia_bio_names_abund) %>%
  mutate(predictor = geo_names[match(parameter, prednames)],
         response = factor(bio_titles_abundance[match(rv, fia_bio_names_abund)], levels = bio_titles_abundance))
spatial_r2s <- spatial_r2s %>%
  mutate(response = factor(c(bio_titles_incidence, bio_titles_abundance)[match(rv, c(fia_bio_names_incid, fia_bio_names_abund))], levels = c(bio_titles_incidence, bio_titles_abundance))) %>%
  left_join(rmse_all)

coefplot_bbs <- coef_bbs_fixed %>%
  group_by(ecoregion) %>%
  do(plot = ggplot(., aes(x = predictor, y = Estimate)) +
    facet_wrap(~ response, scales = 'free_y') +
    geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = 0) +
    geom_text(aes(label = paste('R^2 ==', round(R2, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -1.7, parse = TRUE, data = filter(spatial_r2s, taxon == 'bbs', ecoregion %in% .$ecoregion)) +
    geom_text(aes(label = paste('rRMSE =', round(rRMSE, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -.7, parse = FALSE, data = filter(spatial_r2s, taxon == 'bbs', ecoregion %in% .$ecoregion)) +  
    geom_point() +
    theme_bw() +
    theme(strip.background = element_rect(fill=NA),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle('BBS Geodiversity predictors of Biodiversity at 100 km scale', .$ecoregion[1])
  )

coefplot_fia_incid <- coef_fia_fixed_incid %>%
  group_by(ecoregion) %>%
  do(plot = ggplot(., aes(x = predictor, y = Estimate)) +
       facet_wrap(~ response, scales = 'free_y') +
       geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
       geom_errorbar(aes(ymin = q025, ymax = q975), width = 0) +
       geom_text(aes(label = paste('R^2 ==', round(R2, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -1.7, parse = TRUE, data = filter(spatial_r2s, taxon == 'fia', ecoregion %in% .$ecoregion, !grepl('abundance', response))) +
       geom_text(aes(label = paste('rRMSE =', round(rRMSE, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -.7, parse = FALSE, data = filter(spatial_r2s, taxon == 'fia', ecoregion %in% .$ecoregion, !grepl('abundance', response))) +  
       geom_point() +
       theme_bw() +
       theme(strip.background = element_rect(fill=NA),
             axis.text.x = element_text(angle = 45, hjust = 1)) +
       ggtitle('FIA Geodiversity predictors of Biodiversity at 100 km scale (incidence)', .$ecoregion[1])
  )

coefplot_fia_abund <- coef_fia_fixed_abund %>%
  group_by(ecoregion) %>%
  do(plot = ggplot(., aes(x = predictor, y = Estimate)) +
       facet_wrap(~ response, scales = 'free_y') +
       geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
       geom_errorbar(aes(ymin = q025, ymax = q975), width = 0) +
       geom_text(aes(label = paste('R^2 ==', round(R2, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -1.7, parse = TRUE, data = filter(spatial_r2s, taxon == 'fia', ecoregion %in% .$ecoregion, grepl('abundance', response))) +
       geom_text(aes(label = paste('rRMSE =', round(rRMSE, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -.7, parse = FALSE, data = filter(spatial_r2s, taxon == 'fia', ecoregion %in% .$ecoregion, grepl('abundance', response))) +         geom_point() +
       theme_bw() +
       theme(strip.background = element_rect(fill=NA),
             axis.text.x = element_text(angle = 45, hjust = 1)) +
       ggtitle('FIA Geodiversity predictors of Biodiversity at 100 km scale (abundance)', .$ecoregion[1])
  )

region_names <- c('BCR', 'HUC4', 'TNC')
fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/observed_predicted_plots'

for (i in 1:3) {
  ggsave(file.path(fpfig, paste0('BBS_', region_names[i], '_coefficients.png')), coefplot_bbs$plot[[i]], height = 9, width = 8, dpi = 300)
  ggsave(file.path(fpfig, paste0('FIAincidence_', region_names[i], '_coefficients.png')), coefplot_fia_incid$plot[[i]], height = 9, width = 8, dpi = 300)
  ggsave(file.path(fpfig, paste0('FIAabundance_', region_names[i], '_coefficients.png')), coefplot_fia_abund$plot[[i]], height = 9, width = 8, dpi = 300)
}