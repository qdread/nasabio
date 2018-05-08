# Fixed effect coefficient plots from spatial mixed models
# QDR/NASABioxgeo/27 Apr 2018

# Edit 08 May: do separately for 50 and 100 km, and add the nonspatial models too.
# Edit 03 May: add k-fold RMSE
# Edit 02 May: add RMSE and observed vs predicted plots

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/modelfits' # Local

coef50 <- read.csv(file.path(fp, 'spatial_coef_50k.csv'), stringsAsFactors = FALSE)
coef100 <- read.csv(file.path(fp, 'spatial_coef_100k.csv'), stringsAsFactors = FALSE)
coef50non <- read.csv(file.path(fp, 'nonspatial_coef_50k.csv'), stringsAsFactors = FALSE)
coef100non <- read.csv(file.path(fp, 'nonspatial_coef_100k.csv'), stringsAsFactors = FALSE)

spatial_r2s50 <- read.csv(file.path(fp, 'spatial_r2s_50k.csv'), stringsAsFactors = FALSE)
spatial_r2s100 <- read.csv(file.path(fp, 'spatial_r2s_100k.csv'), stringsAsFactors = FALSE)
nonspatial_r2s50 <- read.csv(file.path(fp, 'nonspatial_r2s_50k.csv'), stringsAsFactors = FALSE)
nonspatial_r2s100 <- read.csv(file.path(fp, 'nonspatial_r2s_100k.csv'), stringsAsFactors = FALSE)

pred50 <- read.csv(file.path(fp, 'spatial_pred_50k.csv'), stringsAsFactors = FALSE)
pred100 <- read.csv(file.path(fp, 'spatial_pred_100k.csv'), stringsAsFactors = FALSE)
pred50non <- read.csv(file.path(fp, 'nonspatial_pred_50k.csv'), stringsAsFactors = FALSE)
pred100non <- read.csv(file.path(fp, 'nonspatial_pred_100k.csv'), stringsAsFactors = FALSE)

kfold_stats50 <- read.csv(file.path(fp, 'spatial_kfold_stats_50k.csv'), stringsAsFactors = FALSE)
kfold_stats100 <- read.csv(file.path(fp, 'spatial_kfold_stats_100k.csv'), stringsAsFactors = FALSE)
kfold_stats50non <- read.csv(file.path(fp, 'nonspatial_kfold_stats_50k.csv'), stringsAsFactors = FALSE)
kfold_stats100non <- read.csv(file.path(fp, 'nonspatial_kfold_stats_100k.csv'), stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)

# Calculation of RMSE
rmse50 <- pred50 %>%
  group_by(taxon, rv, ecoregion) %>%
  summarize(RMSE = sqrt(mean((observed-Estimate)^2)),
            range_obs = diff(range(observed)),
            rRMSE = RMSE/range_obs)
rmse100 <- pred100 %>%
  group_by(taxon, rv, ecoregion) %>%
  summarize(RMSE = sqrt(mean((observed-Estimate)^2)),
            range_obs = diff(range(observed)),
            rRMSE = RMSE/range_obs)
rmse50non <- pred50non %>%
  mutate(ecoregion = 'none') %>%
  group_by(taxon, rv, ecoregion) %>%
  summarize(RMSE = sqrt(mean((observed-Estimate)^2)),
            range_obs = diff(range(observed)),
            rRMSE = RMSE/range_obs)
rmse100non <- pred100non %>%
  mutate(ecoregion = 'none') %>%
  group_by(taxon, rv, ecoregion) %>%
  summarize(RMSE = sqrt(mean((observed-Estimate)^2)),
            range_obs = diff(range(observed)),
            rRMSE = RMSE/range_obs)

rmse_all <- rbind(data.frame(radius = 50, rbind(rmse50, rmse50non)), data.frame(radius = 100, rbind(rmse100, rmse100non)))

prednames50 <- c('elevation_5k_50_sd', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'bio12_5k_50_sd', 'dhi_gpp_5k_50_sd', 'human_footprint_5k_50_mean')
prednames100 <- c('elevation_5k_100_sd', 'bio1_5k_100_mean', 'geological_age_5k_100_diversity', 'soil_type_5k_100_diversity', 'bio12_5k_100_mean', 'bio12_5k_100_sd', 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean')
geo_names <- c('elevation sd','temperature mean','geol age diversity','soil diversity','precip mean','precip sd','gpp sd','footprint mean')
geo_names_order <- c('temperature mean', 'precip mean', 'footprint mean', 'elevation sd', 'precip sd', 'gpp sd', 'geol age diversity', 'soil diversity')

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

# Combine all the data frames for each type of coefficient.

coef_all <- bind_rows(data.frame(radius = 50, coef50), data.frame(radius = 100, coef100),
                      data.frame(ecoregion = 'none', effect = 'fixed', region = NA, coef50non), data.frame(ecoregion = 'none', effect = 'fixed', region = NA, coef100non)) %>%
  mutate(predictor = c(geo_names,geo_names)[match(parameter, c(prednames50, prednames100))],
         response = factor(c(bio_titles,bio_titles_abundance)[match(rv, c(bio_names,fia_bio_names_abund))], levels = c(bio_titles,bio_titles_abundance))) %>%
  mutate(predictor = factor(predictor, levels = geo_names_order),
         ecoregion = factor(ecoregion, levels = c('HUC4','TNC','BCR','none')))

r2_all <- bind_rows(data.frame(radius = 50, spatial_r2s50), data.frame(radius = 100, spatial_r2s100),
                    data.frame(ecoregion = 'none', nonspatial_r2s50), data.frame(ecoregion = 'none', nonspatial_r2s100))

kfold_all <- bind_rows(data.frame(radius = 50, kfold_stats50), data.frame(radius = 100, kfold_stats100),
                    data.frame(ecoregion = 'none', kfold_stats50non), data.frame(ecoregion = 'none', kfold_stats100non))

r2_all <- r2_all %>%
  mutate(response = factor(c(bio_titles_incidence, bio_titles_abundance)[match(rv, c(fia_bio_names_incid, fia_bio_names_abund))], levels = c(bio_titles_incidence, bio_titles_abundance))) %>%
  left_join(rmse_all) %>%
  left_join(kfold_all) %>%
  mutate(kfold_rRMSE = rmse_total/range_obs,
         ecoregion = factor(ecoregion, levels = c('HUC4','TNC','BCR','none')))

coefplot_bbs <- coef_bbs_fixed %>%
  group_by(ecoregion) %>%
  do(plot = ggplot(., aes(x = predictor, y = Estimate)) +
    facet_wrap(~ response, scales = 'free_y') +
    geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
    geom_errorbar(aes(ymin = q025, ymax = q975), width = 0) +
    geom_text(aes(label = paste('R^2 ==', round(R2, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -2.7, parse = TRUE, data = filter(spatial_r2s, taxon == 'bbs', ecoregion %in% .$ecoregion)) +
    geom_text(aes(label = paste('rRMSE =', round(rRMSE, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -1.7, parse = FALSE, data = filter(spatial_r2s, taxon == 'bbs', ecoregion %in% .$ecoregion)) +  
    geom_text(aes(label = paste('out of sample rRMSE =', round(kfold_rRMSE, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -0.7, parse = FALSE, data = filter(spatial_r2s, taxon == 'bbs', ecoregion %in% .$ecoregion)) +    
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
       geom_text(aes(label = paste('R^2 ==', round(R2, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -2.7, parse = TRUE, data = filter(spatial_r2s, taxon == 'fia', ecoregion %in% .$ecoregion, !grepl('abundance', response))) +
       geom_text(aes(label = paste('rRMSE =', round(rRMSE, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -1.7, parse = FALSE, data = filter(spatial_r2s, taxon == 'fia', ecoregion %in% .$ecoregion, !grepl('abundance', response))) +  
       geom_text(aes(label = paste('out of sample rRMSE =', round(kfold_rRMSE, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -.7, parse = FALSE, data = filter(spatial_r2s, taxon == 'fia', ecoregion %in% .$ecoregion, !grepl('abundance', response))) +
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
       geom_text(aes(label = paste('R^2 ==', round(R2, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -2.7, parse = TRUE, data = filter(spatial_r2s, taxon == 'fia', ecoregion %in% .$ecoregion, grepl('abundance', response))) +
       geom_text(aes(label = paste('rRMSE =', round(rRMSE, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -1.7, parse = FALSE, data = filter(spatial_r2s, taxon == 'fia', ecoregion %in% .$ecoregion, grepl('abundance', response))) +       
       geom_text(aes(label = paste('out of sample rRMSE =', round(kfold_rRMSE, 2))), x = -Inf, y = -Inf, hjust = -.1, vjust = -.7, parse = FALSE, data = filter(spatial_r2s, taxon == 'fia', ecoregion %in% .$ecoregion, grepl('abundance', response))) +         
       geom_point() +
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



# All regions on one plot -------------------------------------------------

pd <- position_dodge(width = 0.3)


coefplot_bbs_50 <- coef_all %>%
  filter(taxon == 'bbs', radius == 50, effect == 'fixed', !parameter %in% 'Intercept') %>%
  ggplot(aes(x = predictor, y = Estimate, color = ecoregion)) +
       facet_wrap(~ response, scales = 'free_y') +
       geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
       geom_errorbar(aes(ymin = q025, ymax = q975), width = 0, position = pd) +
       geom_point(position = pd) +
       theme_bw() +
       theme(strip.background = element_rect(fill=NA),
             axis.text.x = element_text(angle = 45, hjust = 1),
             legend.position = 'bottom') +
       scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3','black')) +
       ggtitle('BBS Geodiversity predictors of Biodiversity at 50 km scale')

coefplot_fia_incid_50 <- coef_all %>%
  filter(taxon == 'fia', radius == 50, effect == 'fixed', !parameter %in% 'Intercept', !grepl('abundance',response)) %>%
  ggplot(aes(x = predictor, y = Estimate, color = ecoregion)) +
  facet_wrap(~ response, scales = 'free_y') +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0, position = pd) +
  geom_point(position = pd) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom') +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3','black')) +
  ggtitle('FIA Geodiversity predictors of Biodiversity at 50 km scale', 'incidence-based')

coefplot_fia_abund_50 <- coef_all %>%
  filter(taxon == 'fia', radius == 50, effect == 'fixed', !parameter %in% 'Intercept', grepl('abundance',response)) %>%
  ggplot(aes(x = predictor, y = Estimate, color = ecoregion)) +
  facet_wrap(~ response, scales = 'free_y') +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0, position = pd) +
  geom_point(position = pd) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom') +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3','black')) +
  ggtitle('FIA Geodiversity predictors of Biodiversity at 50 km scale', 'abundance-based')

coefplot_bbs_100 <- coef_all %>%
  filter(taxon == 'bbs', radius == 100, effect == 'fixed', !parameter %in% 'Intercept') %>%
  ggplot(aes(x = predictor, y = Estimate, color = ecoregion)) +
  facet_wrap(~ response, scales = 'free_y') +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0, position = pd) +
  geom_point(position = pd) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom') +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3','black')) +
  ggtitle('BBS Geodiversity predictors of Biodiversity at 100 km scale')

coefplot_fia_incid_100 <- coef_all %>%
  filter(taxon == 'fia', radius == 100, effect == 'fixed', !parameter %in% 'Intercept', !grepl('abundance',response)) %>%
  ggplot(aes(x = predictor, y = Estimate, color = ecoregion)) +
  facet_wrap(~ response, scales = 'free_y') +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0, position = pd) +
  geom_point(position = pd) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom') +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3','black')) +
  ggtitle('FIA Geodiversity predictors of Biodiversity at 100 km scale', 'incidence-based')

coefplot_fia_abund_100 <- coef_all %>%
  filter(taxon == 'fia', radius == 100, effect == 'fixed', !parameter %in% 'Intercept', grepl('abundance',response)) %>%
  ggplot(aes(x = predictor, y = Estimate, color = ecoregion)) +
  facet_wrap(~ response, scales = 'free_y') +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0, position = pd) +
  geom_point(position = pd) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom') +
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3','black')) +
  ggtitle('FIA Geodiversity predictors of Biodiversity at 100 km scale', 'abundance-based')


# Save the plots

fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/observed_predicted_plots'

ggsave(file.path(fpfig, 'BBS_all_coefficients_50km.png'), coefplot_bbs_50, height = 9, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'FIAincidence_all_coefficients_50km.png'), coefplot_fia_incid_50, height = 9, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'FIAabundance_all_coefficients_50km.png'), coefplot_fia_abund_50, height = 9, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'BBS_all_coefficients_100km.png'), coefplot_bbs_100, height = 9, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'FIAincidence_all_coefficients_100km.png'), coefplot_fia_incid_100, height = 9, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'FIAabundance_all_coefficients_100km.png'), coefplot_fia_abund_100, height = 9, width = 8, dpi = 300)


# Just TNC coefficients ---------------------------------------------------

# Add some color to indicate which ones' credible intervals are not zero

coefplot_bbs_tnc50 <- coef_all %>%
  filter(taxon == 'bbs', radius == 50, effect == 'fixed', !parameter %in% 'Intercept', ecoregion == 'TNC') %>%
  mutate(nonzero = q025 > 0 | q975 < 0) %>%
  ggplot(aes(x = predictor, y = Estimate, color = nonzero)) +
  facet_wrap(~ response, scales = 'free_y') +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0) +
  geom_point() +
  scale_color_manual(values = c('black', 'red')) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')

coefplot_fia_incid_tnc50 <- coef_all %>%
  filter(taxon == 'fia', radius == 50, effect == 'fixed', !parameter %in% 'Intercept', !grepl('abundance',response), ecoregion == 'TNC') %>%
  mutate(nonzero = q025 > 0 | q975 < 0) %>%
  ggplot(aes(x = predictor, y = Estimate, color = nonzero)) +
  facet_wrap(~ response, scales = 'free_y') +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(ymin = q025, ymax = q975), width = 0) +
  geom_point() +
  scale_color_manual(values = c('black', 'red')) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none') 

ggsave(file.path(fpfig, 'BBS_TNC_coefficients_50km.png'), coefplot_bbs_tnc50, height = 8, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'FIAincidence_TNC_coefficients_50km.png'), coefplot_fia_incid_tnc50, height = 8, width = 8, dpi = 300)


# Plot showing RMSEs --------------------------------------------------------

rmseplot_bbs_50 <- r2_all %>% 
  filter(taxon == 'bbs', radius == 50) %>%
  ggplot(aes(x = ecoregion)) +
  geom_point(aes(y = rRMSE)) +
  geom_point(aes(y = kfold_rRMSE), color = 'red') +
  facet_wrap(~ response) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom') +
  ggtitle('BBS model performance, 50 km scale')

rmseplot_bbs_100 <- r2_all %>% 
  filter(taxon == 'bbs', radius == 100) %>%
  ggplot(aes(x = ecoregion)) +
  geom_point(aes(y = rRMSE)) +
  geom_point(aes(y = kfold_rRMSE), color = 'red') +
  facet_wrap(~ response) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom') +
  ggtitle('BBS model performance, 100 km scale')

rmseplot_fia_incid_50 <- r2_all %>% 
  filter(taxon == 'fia', radius == 50, !grepl('abundance', response)) %>%
  ggplot(aes(x = ecoregion)) +
  geom_point(aes(y = rRMSE)) +
  geom_point(aes(y = kfold_rRMSE), color = 'red') +
  facet_wrap(~ response) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom') +
  ggtitle('FIA model performance, 50 km scale', 'incidence-based')

rmseplot_fia_incid_100 <- r2_all %>% 
  filter(taxon == 'fia', radius == 100, !grepl('abundance', response)) %>%
  ggplot(aes(x = ecoregion)) +
  geom_point(aes(y = rRMSE)) +
  geom_point(aes(y = kfold_rRMSE), color = 'red') +
  facet_wrap(~ response) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom') +
  ggtitle('FIA model performance, 100 km scale', 'incidence-based')

rmseplot_fia_abund_50 <- r2_all %>% 
  filter(taxon == 'fia', radius == 50, grepl('abundance', response)) %>%
  ggplot(aes(x = ecoregion)) +
  geom_point(aes(y = rRMSE)) +
  geom_point(aes(y = kfold_rRMSE), color = 'red') +
  facet_wrap(~ response) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom') +
  ggtitle('FIA model performance, 50 km scale', 'abundance-based')

rmseplot_fia_abund_100 <- r2_all %>% 
  filter(taxon == 'fia', radius == 100, grepl('abundance', response)) %>%
  ggplot(aes(x = ecoregion)) +
  geom_point(aes(y = rRMSE)) +
  geom_point(aes(y = kfold_rRMSE), color = 'red') +
  facet_wrap(~ response) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'bottom') +
  ggtitle('FIA model performance, 100 km scale', 'abundance-based')

# rmse for only TNC ecoregions for BBS
rmseplot_bbs_tnc_50 <- r2_all %>% 
  filter(taxon == 'bbs', radius == 50, ecoregion == 'TNC') %>%
  ggplot(aes(x = response)) +
  geom_point(aes(y = rRMSE)) +
  geom_point(aes(y = kfold_rRMSE), color = 'red') +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.17), expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('BBS model performance, 50 km scale')

rmseplot_fia_tnc_50 <- r2_all %>% 
  filter(taxon == 'fia', radius == 50, ecoregion == 'TNC', !grepl('abundance', response)) %>%
  ggplot(aes(x = response)) +
  geom_point(aes(y = rRMSE)) +
  geom_point(aes(y = kfold_rRMSE), color = 'red') +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.17), expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('FIA model performance, 50 km scale')

# rmseplot for both
rmseplot_both_tnc_50 <- r2_all %>% 
  filter(radius == 50, ecoregion == 'TNC', !grepl('abundance', response)) %>%
  ggplot(aes(x = response)) +
  facet_grid(. ~ taxon, labeller = labeller(taxon = c(bbs = 'birds', fia = 'trees'))) +
  geom_point(aes(y = rRMSE)) +
  geom_point(aes(y = kfold_rRMSE), color = 'red') +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.17), expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = NA)) 

# Save the plots

fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/observed_predicted_plots'

ggsave(file.path(fpfig, 'BBS_performance_50km.png'), rmseplot_bbs_50, height = 9, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'FIAincidence_performance_50km.png'), rmseplot_fia_incid_50, height = 9, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'FIAabundance_performance_50km.png'), rmseplot_fia_abund_50, height = 9, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'BBS_performance_100km.png'), rmseplot_bbs_100, height = 9, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'FIAincidence_performance_100km.png'), rmseplot_fia_incid_100, height = 9, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'FIAabundance_performance_100km.png'), rmseplot_fia_abund_100, height = 9, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'BBS_performance_TNC.png'), rmseplot_bbs_tnc_50, height = 5, width = 5, dpi = 300)
ggsave(file.path(fpfig, 'FIA_performance_TNC.png'), rmseplot_fia_tnc_50, height = 5, width = 5, dpi = 300)
ggsave(file.path(fpfig, 'both_performance_TNC.png'), rmseplot_both_tnc_50, height = 4, width = 7, dpi = 300)
