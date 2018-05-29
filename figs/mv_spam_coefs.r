# Fixed effect coefficient plots from spatial mixed models (multivariate version)
# QDR/NASABioxgeo/28 May 2018


# Load and combine data ---------------------------------------------------


fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/modelfits' # Local

model_coef <- read.csv(file.path(fp, 'multivariate_spatial_coef.csv'), stringsAsFactors = FALSE)
model_pred <- read.csv(file.path(fp, 'multivariate_spatial_pred.csv'), stringsAsFactors = FALSE)
model_rmse <- read.csv(file.path(fp, 'multivariate_spatial_rmse.csv'), stringsAsFactors = FALSE)
kfold_pred <- read.csv(file.path(fp, 'multivariate_kfold_pred.csv'), stringsAsFactors = FALSE)
kfold_rmse <- read.csv(file.path(fp, 'multivariate_kfold_rmse.csv'), stringsAsFactors = FALSE)


library(dplyr)
library(ggplot2)
library(reshape2)

prednames50 <- c('elevation_5k_50_sd', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'bio12_5k_50_sd', 'dhi_gpp_5k_50_sd')
geo_names <- c('elevation sd','temperature mean','geol. age diversity','soil diversity','precip. mean','precip. sd','GPP sd')
geo_names_order <- c('temperature mean', 'precip. mean', 'elevation sd', 'precip. sd', 'GPP sd', 'geol. age diversity', 'soil diversity')

bio_titles <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'beta PD', 'gamma PD', 'alpha FD', 'beta FD', 'gamma FD')
bio_names <- c("alpha_richness", "beta_td_sorensen_pa", "gamma_richness",
               "alpha_phy_pa", "beta_phy_pa", "gamma_phy_pa", 
               "alpha_func_pa", "beta_func_pa", "gamma_func_pa")

raw_bio_names <- cbind(expand.grid(rv=c('alpha','beta','gamma'),variable=c('rmse_y1','rmse_y2','rmse_y3')),
                       name = bio_names)

# K-fold RMSE is still separated by fold, so compute the total RMSE for each variable.
# Combine full-model and k-fold RMSEs.
all_rmse <- kfold_rmse %>%
  select(-kfoldic, -kfoldic_se) %>%
  melt(id.vars = c('taxon','rv','ecoregion','fold'), value.name = 'kfold_RMSE') %>%
  mutate(variable = raw_bio_names$name[match(paste(rv, variable), paste(raw_bio_names$rv, raw_bio_names$variable))]) %>%
  group_by(taxon, rv, ecoregion, variable) %>%
  summarize(kfold_RMSE = sqrt(mean(kfold_RMSE^2))) %>%
  rename(response = variable) %>%
  right_join(model_rmse) %>%
  mutate(kfold_rRMSE = kfold_RMSE/range_obs,
         response = factor(bio_titles[match(response, bio_names)], levels = bio_titles))

# Reshape coefficient plot and relabel it
all_coef <- model_coef %>%
  filter(effect == 'fixed', !parameter %in% 'Intercept') %>%
  dcast(taxon + rv + response + parameter ~ stat) %>%
  mutate(predictor = factor(geo_names[match(parameter, prednames50)], levels = geo_names_order),
         response = factor(bio_titles[match(response, bio_names)], levels = bio_titles)) 
  

# Coefficient plots -------------------------------------------------

fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/observed_predicted_plots'

# Add some color to indicate which ones' credible intervals are not zero

coefplot_bbs <- all_coef %>%
  filter(taxon == 'bbs') %>%
  mutate(nonzero = Q2.5 > 0 | Q97.5 < 0) %>%
  ggplot(aes(x = predictor, y = Estimate, color = nonzero)) +
  facet_wrap(~ response, scales = 'free_y') +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0) +
  geom_point() +
  scale_color_manual(values = c('black', 'red')) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')

coefplot_fia <- all_coef %>%
  filter(taxon == 'fia') %>%
  mutate(nonzero = Q2.5 > 0 | Q97.5 < 0) %>%
  ggplot(aes(x = predictor, y = Estimate, color = nonzero)) +
  facet_wrap(~ response, scales = 'free_y') +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(ymin = Q2.5, ymax = Q97.5), width = 0) +
  geom_point() +
  scale_color_manual(values = c('black', 'red')) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')


ggsave(file.path(fpfig, 'BBS_multivariate_coef.png'), coefplot_bbs, height = 8, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'FIA_multivariate_coef.png'), coefplot_fia, height = 8, width = 8, dpi = 300)


# Plot showing RMSEs --------------------------------------------------------

# rmse for only TNC ecoregions for BBS
rmseplot_bbs <- all_rmse %>% 
  filter(taxon == 'bbs') %>%
  ggplot(aes(x = response)) +
  geom_point(aes(y = rRMSE)) +
  geom_point(aes(y = kfold_rRMSE), color = 'red') +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.17), expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('BBS model performance')

rmseplot_fia <- all_rmse %>% 
  filter(taxon == 'fia') %>%
  ggplot(aes(x = response)) +
  geom_point(aes(y = rRMSE)) +
  geom_point(aes(y = kfold_rRMSE), color = 'red') +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.17), expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('FIA model performance')

# rmseplot for both
rmseplot_both <- all_rmse %>% 
  ggplot(aes(x = response)) +
  facet_grid(. ~ taxon, labeller = labeller(taxon = c(bbs = 'birds', fia = 'trees'))) +
  geom_point(aes(y = rRMSE)) +
  geom_point(aes(y = kfold_rRMSE), color = 'red') +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.17), expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = NA)) 

# Save the plots

ggsave(file.path(fpfig, 'BBS_performance_multivariate.png'), rmseplot_bbs, height = 5, width = 5, dpi = 300)
ggsave(file.path(fpfig, 'FIA_performance_multivariate.png'), rmseplot_fia, height = 5, width = 5, dpi = 300)
ggsave(file.path(fpfig, 'both_performance_multivariate.png'), rmseplot_both, height = 4, width = 7, dpi = 300)


