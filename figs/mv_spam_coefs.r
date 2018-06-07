# Fixed effect coefficient plots from spatial mixed models (multivariate version)
# QDR/NASABioxgeo/28 May 2018


# Load and combine data ---------------------------------------------------


fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/modelfits' # Local

model_coef <- read.csv(file.path(fp, 'multivariate_spatial_coef.csv'), stringsAsFactors = FALSE)
model_pred <- read.csv(file.path(fp, 'multivariate_spatial_pred.csv'), stringsAsFactors = FALSE)
model_rmse <- read.csv(file.path(fp, 'multivariate_spatial_rmse.csv'), stringsAsFactors = FALSE)
model_r2 <- read.csv(file.path(fp, 'multivariate_spatial_r2.csv'), stringsAsFactors = FALSE)
model_coef_var <- read.csv(file.path(fp, 'multivariate_spatial_coef_variation.csv'), stringsAsFactors = FALSE)
#kfold_pred <- read.csv(file.path(fp, 'multivariate_kfold_pred.csv'), stringsAsFactors = FALSE)
kfold_rmse <- read.csv(file.path(fp, 'multivariate_kfold_rmse.csv'), stringsAsFactors = FALSE)


library(dplyr)
library(ggplot2)
library(reshape2)
library(purrr)

prednames50 <- c('elevation_5k_50_sd', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'bio12_5k_50_sd', 'dhi_gpp_5k_50_sd')
geo_names <- c('elevation sd','temperature mean','geol. age diversity','soil diversity','precip. mean','precip. sd','GPP sd')
geo_names_order <- c('temperature mean', 'precip. mean', 'elevation sd', 'precip. sd', 'GPP sd', 'geol. age diversity', 'soil diversity')

bio_titles <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'beta PD', 'gamma PD', 'alpha FD', 'beta FD', 'gamma FD')
bio_names <- c("alpha_richness", "beta_td_sorensen_pa", "gamma_richness",
               "alpha_phy_pa", "beta_phy_pa", "gamma_phy_pa", 
               "alpha_func_pa", "beta_func_pa", "gamma_func_pa")

raw_bio_names <- cbind(expand.grid(rv=c('alpha','beta','gamma'),variable=c('rmse_y1','rmse_y2','rmse_y3')),
                       name = bio_names)

# Combine full-model and k-fold RMSEs.
# Include RMSE from each fold so we can see variability due to folds.

all_rmse <- kfold_rmse %>%
  mutate(response = bio_names[match(response, gsub('_','',bio_names))]) %>%
  filter(is.na(fold)) %>%
  select(-kfoldic, -kfoldic_se, -fold) %>%
  right_join(model_rmse) %>%
  left_join(model_r2 %>% rename(r2 = Estimate, r2_error = Est.Error, r2_q025 = Q2.5, r2_q975 = Q97.5)) %>%
  mutate_at(vars(starts_with('kfold_RMSE')), funs(relative = ./range_obs)) %>%
  mutate(response = factor(bio_titles[match(response, bio_names)], levels = bio_titles))

# Reshape coefficient plot and relabel it
all_coef <- model_coef %>%
  filter(effect == 'fixed', !parameter %in% 'Intercept') %>%
  dcast(taxon + rv + response + parameter ~ stat) %>%
  mutate(predictor = factor(geo_names[match(parameter, prednames50)], levels = geo_names_order),
         response = factor(bio_titles[match(response, bio_names)], levels = bio_titles),
         flavor = map_chr(strsplit(as.character(response), ' '), 2) %>%
           factor(levels = c('TD','PD','FD'), labels = c('taxonomic', 'phylogenetic', 'functional')))

# Relabel data frame of spatial variability metrics  
model_coef_var <- model_coef_var %>%
  mutate(predictor = factor(geo_names[match(parameter, prednames50)], levels = geo_names_order),
         response = factor(bio_titles[match(response, bio_names)], levels = bio_titles),
         flavor = map_chr(strsplit(as.character(response), ' '), 2) %>%
           factor(levels = c('TD','PD','FD'), labels = c('taxonomic', 'phylogenetic', 'functional')))

# Coefficient plots -------------------------------------------------

fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/multivariate_maps_figs'

# Add some color to indicate which ones' credible intervals are not zero
# Also shade the climate mean region with a gray rectangle

coefdat_bbs <- all_coef %>%
  filter(taxon == 'bbs') %>%
  mutate(nonzero = Q2.5 > 0 | Q97.5 < 0)
coefplot_bbs <- ggplot(coefdat_bbs) +
  geom_rect(xmin=0, xmax=2.5, ymin=-Inf, ymax=Inf, fill = 'gray90') +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(x = predictor, ymin = Q2.5, ymax = Q97.5, color = nonzero), width = 0) +
  geom_point(aes(x = predictor, y = Estimate, color = nonzero)) +
  facet_grid(rv ~ flavor) +
  scale_y_continuous(name = 'coefficient estimate', limits = c(-0.73, 0.73), expand = c(0,0)) +
  scale_color_manual(values = c('black', 'red')) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')

coefdat_fia <- all_coef %>%
  filter(taxon == 'fia') %>%
  mutate(nonzero = Q2.5 > 0 | Q97.5 < 0)
coefplot_fia <- ggplot(coefdat_fia) +
  geom_rect(xmin=0, xmax=2.5, ymin=-Inf, ymax=Inf, fill = 'gray90') +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(x = predictor, ymin = Q2.5, ymax = Q97.5, color = nonzero), width = 0) +
  geom_point(aes(x = predictor, y = Estimate, color = nonzero)) +
  facet_grid(rv ~ flavor) +
  scale_color_manual(values = c('black', 'red')) +
  scale_y_continuous(name = 'coefficient estimate', limits = c(-0.73, 0.73), expand = c(0,0)) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')


ggsave(file.path(fpfig, 'BBS_multivariate_coef.png'), coefplot_bbs, height = 8, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'FIA_multivariate_coef.png'), coefplot_fia, height = 8, width = 8, dpi = 300)

# Edit 6 June: sideways coefficient plot.
coef_bbs_sideways <- coefplot_bbs + 
  coord_flip() +
  theme(axis.text.x = element_text(angle=0, hjust=0.5)) 
coef_fia_sideways <- coefplot_fia + 
  coord_flip() +
  theme(axis.text.x = element_text(angle=0, hjust=0.5)) 

ggsave(file.path(fpfig, 'BBS_multivariate_coef_sideways.png'), coef_bbs_sideways, height = 8, width = 8, dpi = 300)
ggsave(file.path(fpfig, 'FIA_multivariate_coef_sideways.png'), coef_fia_sideways, height = 8, width = 8, dpi = 300)


# Plot of spatial variability ---------------------------------------------

pd = position_dodge(width = 0.5)
coefvar_plot <- ggplot(model_coef_var %>% 
                         filter(!is.na(predictor)) %>%
                         mutate(taxon = factor(taxon,labels=c('birds','trees')))) +
  geom_rect(xmin=0, xmax=2.5, ymin=-Inf, ymax=Inf, fill = 'gray90') +
  geom_col(aes(x = predictor, y = cv, fill = taxon, group = taxon), position = pd, width = 0.5) +
  facet_grid(rv ~ flavor) +
  scale_fill_manual(values = c('blue', 'skyblue')) +
  scale_y_continuous(name = 'spatial coefficient of variation', limits = c(0, 15), expand = c(0,0), breaks=c(0,5,10)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(0.91, 0.9))

coefvar_sideways <- coefvar_plot + 
  coord_flip() +
  theme(axis.text.x = element_text(angle=0, hjust=0.5),
        legend.background = element_rect(color = 'black')) 
ggsave(file.path(fpfig, 'both_multivariate_coefvar_sideways.png'), coefvar_sideways, height = 8, width = 8, dpi = 300)

# Plot showing RMSEs --------------------------------------------------------

pn1 <- position_nudge(x = -0.06, y = 0)
pn2 <- position_nudge(x = 0.06, y = 0)
# rmseplot for both
rmseplot_both <- all_rmse %>% 
  ggplot(aes(x = response)) +
  facet_grid(. ~ taxon, labeller = labeller(taxon = c(bbs = 'birds', fia = 'trees'))) +
  geom_errorbar(aes(ymin = RMSE_q025, ymax = RMSE_q975, y = RMSE_mean), width = 0, position = pn1) +
  geom_errorbar(aes(ymin = kfold_RMSE_q025, ymax = kfold_RMSE_q975, y = RMSE_mean), width = 0, color = 'red', position = pn2) +
  geom_point(aes(y = RMSE_mean), position = pn1) +
  geom_point(aes(y = kfold_RMSE_mean), color = 'red', position = pn2) +
  geom_text(aes(label = round(r2, 2)), y = -Inf, vjust = -0.2, fontface = 3, size = 3) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1.3), expand = c(0,0), name = 'root mean squared error') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = NA)) 

# Save the plots

ggsave(file.path(fpfig, 'both_performance_multivariate.png'), rmseplot_both, height = 4, width = 7, dpi = 300)


