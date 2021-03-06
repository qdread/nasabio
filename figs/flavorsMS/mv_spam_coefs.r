# Fixed effect coefficient plots from spatial mixed models (multivariate version)
# QDR/NASABioxgeo/28 May 2018

### THESE ARE THE ACTUAL PUBLICATION FIGURES FOR THE "FLAVORS" MANUSCRIPT!

# Modified 24 Dec 2019: make figures legible in grayscale.
# Modified 20 Dec 2019: also output PDFs
# Modified 08 May 2019: add WAIC to summary stats table
# Modified 08 Jan 2019: change the blue color scheme and fix the mismatch in the spatial var color scheme.
# Modified 18 June: Add the new null models.
# Modified 15 November: Make a few visual changes to figs 5 and 6 for better readability.

# Load and combine data ---------------------------------------------------


fp <- '~/Dropbox/Q/projects/nasabiodiv/modelfits' # Local

model_coef <- read.csv(file.path(fp, 'multivariate_spatial_coef.csv'), stringsAsFactors = FALSE)
model_pred <- read.csv(file.path(fp, 'multivariate_spatial_pred.csv'), stringsAsFactors = FALSE)
model_rmse <- read.csv(file.path(fp, 'multivariate_spatial_rmse.csv'), stringsAsFactors = FALSE)
model_r2 <- read.csv(file.path(fp, 'multivariate_spatial_r2.csv'), stringsAsFactors = FALSE)
model_coef_var <- read.csv(file.path(fp, 'multivariate_spatial_coef_variation_corrected.csv'), stringsAsFactors = FALSE)
model_waic <- read.csv(file.path(fp, 'multivariate_spatial_waic.csv'), stringsAsFactors = FALSE)
kfold_rmse <- read.csv(file.path(fp, 'multivariate_kfold_rmse.csv'), stringsAsFactors = FALSE)


library(dplyr)
library(ggplot2)
library(reshape2)
library(purrr)

prednames50 <- c('elevation_5k_tri_50_mean', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'dhi_gpp_5k_tri_50_mean')
geo_names <- c('elevation diversity','temperature mean','geol. age diversity','soil diversity','precip. mean','GPP diversity')
geo_names_order <- c('temperature mean', 'precip. mean', 'elevation diversity', 'GPP diversity', 'geol. age diversity', 'soil diversity')

bio_titles <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'beta PD', 'gamma PD', 'alpha FD', 'beta FD', 'gamma FD')
bio_names <- c("alpha_richness", "beta_td_sorensen_pa", "gamma_richness",
               "alpha_phy_pa", "beta_phy_pa", "gamma_phy_pa", 
               "alpha_func_pa", "beta_func_pa", "gamma_func_pa")

twocolors <- c('black', 'gray60')

# Combine full-model and k-fold RMSEs.
# Include RMSE from each fold so we can see variability due to folds.
# Edited 02 May 2019: correct the data wrangling code because the output of the new k-fold is slightly different.

all_rmse <- kfold_rmse %>%
  rename_if(is.numeric, ~ paste0('kfold_', .x)) %>%
  right_join(model_rmse) %>%
  left_join(model_r2 %>% select(-fold) %>% rename(r2 = Estimate, r2_error = Est.Error, r2_q025 = Q2.5, r2_q975 = Q97.5)) %>%
  mutate(response = factor(bio_titles[match(response, bio_names)], levels = bio_titles),
         flavor = map_chr(strsplit(as.character(response), ' '), 2) %>%
           factor(levels = c('TD','PD','FD'), labels = c('taxonomic', 'phylogenetic', 'functional')))
  

# Reshape coefficient plot and relabel it
all_coef <- model_coef %>%
  filter(effect == 'fixed', !parameter %in% 'Intercept') %>%
  dcast(taxon + rv + model + response + parameter ~ stat) %>%
  mutate(predictor = factor(geo_names[match(parameter, prednames50)], levels = geo_names_order),
         response = factor(bio_titles[match(response, bio_names)], levels = bio_titles),
         flavor = map_chr(strsplit(as.character(response), ' '), 2) %>%
           factor(levels = c('TD','PD','FD'), labels = c('taxonomic', 'phylogenetic', 'functional')))

# Relabel data frame of spatial variability metrics  
model_coef_var <- model_coef_var %>%
  mutate(predictor = factor(geo_names[match(parameter, prednames50)], levels = geo_names_order),
         response = factor(bio_titles[match(response, gsub('_', '', bio_names))], levels = bio_titles),
         flavor = map_chr(strsplit(as.character(response), ' '), 2) %>%
           factor(levels = c('TD','PD','FD'), labels = c('taxonomic', 'phylogenetic', 'functional'))) %>%
  rename(coef_var = Estimate)

# Coefficient plots -------------------------------------------------

fpfig <- '~/google_drive/NASABiodiversityWG/Figures/multivariate_maps_figs/07may2019' # updated file path to put new figs in a separate location

# Add some color to indicate which ones' credible intervals are not zero
# Also shade the climate mean region with a gray rectangle

coefdat_bbs <- all_coef %>%
  filter(taxon == 'bbs') %>%
  mutate(nonzero = Q2.5 > 0 | Q97.5 < 0)
coefplot_bbs <- ggplot(coefdat_bbs %>% filter(model=='full')) +
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
coefplot_fia <- ggplot(coefdat_fia %>% filter(model=='full')) +
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

# Edit 25 July 2019: change legend position, add scale shape, resize so things are bigger
# Edit 18 June: plot of coefficients with both, to compare.
# THIS IS FIG 5 IN THE MANUSCRIPT.
coefdat_both <- all_coef %>%
  filter(model == 'full') %>%
  mutate(nonzero = Q2.5 > 0 | Q97.5 < 0)
pd <- position_dodge(width = 0.12)
coefplot_both <- ggplot(coefdat_both %>% mutate(taxon = factor(taxon, labels = c('birds','trees')))) +
  geom_rect(xmin=0, xmax=2.5, ymin=-Inf, ymax=Inf, fill = 'gray95') +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(x = predictor, ymin = Q2.5, ymax = Q97.5, color = taxon, group = taxon, linetype = nonzero), width = 0, position=pd) +
  geom_point(aes(x = predictor, y = Estimate, size = nonzero, color = taxon, fill = nonzero, group = taxon, shape = taxon), position=pd) +
  facet_grid(rv ~ flavor) +
  scale_color_manual(values = twocolors) +
  scale_fill_manual(values = c('black', 'indianred1')) +
  scale_size_manual(values = c(1.5, 2)) +
  scale_shape_manual(values = c(24, 21)) +
  scale_linetype_manual(values = c(3, 1)) +
  scale_y_continuous(name = 'coefficient estimate', limits = c(-0.73, 0.73), expand = c(0,0)) +
  theme_bw() +
  theme(strip.background = element_rect(fill=NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')

coef_both_sideways <- coefplot_both + 
  coord_flip() +
  theme(axis.text.x = element_text(angle=0, hjust=0.5), 
        legend.background = element_rect(color = 'black'),
        legend.position = c(0.92,0.93)) +
  guides(shape = guide_legend(override.aes = list(fill = twocolors)), fill = 'none', size = 'none', linetype = 'none')



ggsave('~/google_drive/NASABiodiversityWG/FlavorsOfDiversityPaper/figures_pdf/fig5.pdf', coef_both_sideways, height = 8, width = 8)
ggsave(file.path(fpfig, 'both_multivariate_coef_sideways.png'), coef_both_sideways, height = 8, width = 8, dpi = 300)


# Plot of spatial variability ---------------------------------------------

pd = position_dodge(width = 0.5)
coefvar_plot <- ggplot(model_coef_var %>% 
                         filter(!is.na(predictor)) %>%
                         mutate(taxon = factor(taxon,levels=c('fia','bbs'),labels=c('trees','birds')),
                                response = map_chr(strsplit(as.character(response), ' '), 1))) +
  geom_rect(xmin=0, xmax=2.5, ymin=-Inf, ymax=Inf, fill = 'gray90') +
  geom_col(aes(x = predictor, y = coef_var, fill = taxon, group = taxon), position = pd, width = 0.5) +
  geom_errorbar(aes(x = predictor, ymin = q025, ymax = q975, group = taxon), position = pd, width = 0.15, color = 'gray20') +
  facet_grid(response ~ flavor) +
  scale_fill_manual(values = rev(twocolors)) +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_y_continuous(name = 'spatial variability of relationship', limits = c(0, 1.16), expand = c(0,0)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(0.92, 0.92))

coefvar_sideways <- coefvar_plot + 
  coord_flip() +
  theme(axis.text.x = element_text(angle=0, hjust=0.5),
        legend.background = element_rect(color = 'black')) 

ggsave(file.path(fpfig, 'both_multivariate_coefvar_sideways.png'), coefvar_sideways, height = 8, width = 8, dpi = 300)

# Plot showing RMSEs --------------------------------------------------------

pn1 <- position_nudge(x = -0.06, y = 0)
pn2 <- position_nudge(x = 0.06, y = 0)
# rmseplot for both
# Comparison of RMSE and R-squared among models.
rmseplot_both <- all_rmse %>% 
  filter(model == 'full') %>%
  ggplot(aes(x = response)) +
  facet_grid(. ~ taxon, labeller = labeller(taxon = c(bbs = 'birds', fia = 'trees'))) +
  geom_errorbar(aes(ymin = RMSE_q025_relative, ymax = RMSE_q975_relative), width = 0, position = pn1) +
  geom_errorbar(aes(ymin = kfold_RMSE_q025_relative, ymax = kfold_RMSE_q975_relative), width = 0, color = 'red', position = pn2) +
  geom_point(aes(y = RMSE_mean_relative), position = pn1) +
  geom_point(aes(y = kfold_RMSE_mean_relative), color = 'red', position = pn2) +
  geom_text(aes(label = round(r2, 2)), y = -Inf, vjust = -0.2, fontface = 3, size = 3) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.65), expand = c(0,0), name = 'relative root mean squared error') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = NA)) 

# Save the plots

ggsave(file.path(fpfig, 'both_performance_multivariate.png'), rmseplot_both, height = 4, width = 7, dpi = 300)

# Edit 18 June: plot comparing RMSEs and R-squared for the 3 model types
# Edit 08 Aug: add geodiversity-only to this

all_rmse <- all_rmse %>%
  mutate(model = factor(model, levels=c('space','climate','geo','full'), labels=c('space only', 'space+climate','space+geodiversity','space+climate+geodiversity')))

pd <- position_dodge(width = 0.05)
rmseplot_bymodel_bird <- all_rmse %>% 
  filter(taxon == 'bbs') %>%
  ggplot(aes(x = rv, color = model, group = model)) +
  facet_grid(. ~ flavor, switch = 'x') +
  geom_errorbar(aes(ymin = RMSE_q025_relative, ymax = RMSE_q975_relative), width = 0, position = pd) +
  geom_point(aes(y = RMSE_mean_relative), position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0,0), name = 'relative root mean squared error') +
  scale_x_discrete(name = 'response') +
  ggtitle('birds') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(0, 'lines'),
        legend.position = 'none')
rmseplot_bymodel_tree <- all_rmse %>% 
  filter(taxon == 'fia') %>%
  ggplot(aes(x = rv, color = model, group = model)) +
  facet_grid(. ~ flavor, switch = 'x') +
  geom_errorbar(aes(ymin = RMSE_q025_relative, ymax = RMSE_q975_relative), width = 0, position = pd) +
  geom_point(aes(y = RMSE_mean_relative), position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.32), expand = c(0,0), name = 'relative root mean squared error') +
  scale_x_discrete(name = 'response') +
  ggtitle('trees') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(0, 'lines'),
        legend.position = c(0.5, 0.2),
        legend.background = element_rect(color = 'black'),
        legend.text = element_text(size = 6.5))


kfold_rmseplot_bymodel_bird <- all_rmse %>% 
  filter(taxon == 'bbs') %>%
  ggplot(aes(x = rv, color = model, group = model)) +
  facet_grid(. ~ flavor, switch = 'x') +
  geom_errorbar(aes(ymin = kfold_RMSE_q025_relative, ymax = kfold_RMSE_q975_relative), width = 0, position = pd) +
  geom_point(aes(y = kfold_RMSE_mean_relative), position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.7), expand = c(0,0), name = 'CV relative root mean squared error') +
  scale_x_discrete(name = 'response') +
  ggtitle('birds') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(0, 'lines'),
        legend.position = 'none')
kfold_rmseplot_bymodel_tree <- all_rmse %>% 
  filter(taxon == 'fia') %>%
  ggplot(aes(x = rv, color = model, group = model)) +
  facet_grid(. ~ flavor, switch = 'x') +
  geom_errorbar(aes(ymin = kfold_RMSE_q025_relative, ymax = kfold_RMSE_q975_relative), width = 0, position = pd) +
  geom_point(aes(y = kfold_RMSE_mean_relative), position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.7), expand = c(0,0), name = 'CV relative root mean squared error') +
  scale_x_discrete(name = 'response') +
  ggtitle('trees') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(0, 'lines'),
        legend.position = c(0.5, 0.2),
        legend.background = element_rect(color = 'black'),
        legend.text = element_text(size = 6.5))

# RMSE plot as 2-way facet
# Update 20 May 2019: add R-squared.
# This is fig 6 in the REVISED manuscript.
# Edit 25 July 2019: include higher jitter and make points bigger
pd <- position_dodge(width = 0.15)

rmseplot_bymodel_2wayfacet <- all_rmse %>% 
  ggplot(aes(x = rv, color = model, group = model)) +
  facet_grid(taxon ~ flavor, switch = 'x', labeller = labeller(taxon = c(bbs = 'birds', fia = 'trees'))) +
  geom_errorbar(aes(ymin = RMSE_q025_relative, ymax = RMSE_q975_relative), width = 0, position = pd) +
  geom_point(aes(y = RMSE_mean_relative), position = pd, size = 2) +
  geom_text(aes(y = -Inf, label = round(r2, 2)), data = all_rmse %>% filter(model == 'space+climate+geodiversity'), color = 'black', fontface = 'italic', vjust = -0.2) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.315), name = 'relative root mean squared error', expand = c(0,0)) +
  scale_x_discrete(name = 'response') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(0, 'lines'),
        legend.position = 'bottom',
        legend.background = element_rect(color = 'black'),
        legend.text = element_text(size = 8))

kfold_rmseplot_bymodel_2wayfacet <- all_rmse %>% 
  ggplot(aes(x = rv, color = model, group = model, shape = model)) +
  facet_grid(taxon ~ flavor, switch = 'x', labeller = labeller(taxon = c(bbs = 'birds', fia = 'trees'))) +
  geom_errorbar(aes(ymin = kfold_RMSE_q025_relative, ymax = kfold_RMSE_q975_relative), width = 0, position = pd) +
  geom_point(aes(y = kfold_RMSE_mean_relative), position = pd, size = 2) +
  theme_bw() +
  scale_shape_manual(values = 21:24) +
  scale_y_continuous(limits = c(0, 0.7), name = 'CV relative root mean squared error', expand = c(0,0)) +
  scale_x_discrete(name = 'response') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(0, 'lines'),
        legend.position = 'bottom',
        legend.background = element_rect(color = 'black'),
        legend.text = element_text(size = 8),
        panel.grid.major.x = element_blank())

# Just compare the predictive power of the full model, putting birds and trees on the same one.
# THIS IS FIG 6 IN THE MANUSCRIPT (or it was in the old version)
rmseplot_taxacolor <- all_rmse %>% 
  filter(model == 'space+climate+geodiversity') %>%
  mutate(taxon = factor(taxon,labels=c('birds','trees'))) %>%
  ggplot(aes(x = rv)) +
  facet_grid(. ~ flavor, switch = 'x') +
  geom_errorbar(aes(ymin = RMSE_q025_relative, ymax = RMSE_q975_relative, color = taxon, group = taxon), width = 0, position = pd) +
  geom_point(aes(y = RMSE_mean_relative, color = taxon, group = taxon), position = pd, size = 2) +
  theme_bw() +
  scale_color_manual(values = twocolors) +
  scale_y_continuous(limits = c(0, 0.31), expand = c(0,0), name = 'relative root mean squared error') +
  scale_x_discrete(name = 'response') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(0, 'lines'),
        legend.position = c(0.91, 0.8),
        legend.background = element_rect(color = 'black'))

# Save plots
library(gridExtra)
png(file.path(fpfig, 'both_fittedrmse_allmodels.png'), height = 8, width = 7, res = 400, units = 'in')
  grid.arrange(rmseplot_bymodel_bird, rmseplot_bymodel_tree, nrow = 2)
dev.off()

png(file.path(fpfig, 'both_lolormse_allmodels.png'), height = 8, width = 7, res = 400, units = 'in')
  grid.arrange(kfold_rmseplot_bymodel_bird, kfold_rmseplot_bymodel_tree, nrow = 2)
dev.off()

ggsave(file.path(fpfig, 'both_fittedrmse_fullonly.png'), rmseplot_taxacolor, height = 4, width = 7, dpi = 400)
ggsave(file.path(fpfig, 'both_fittedrmse_allmodels_2wayfacet.png'), rmseplot_bymodel_2wayfacet, height = 6, width = 7, dpi = 400)
ggsave(file.path(fpfig, 'both_lolormse_allmodels_2wayfacet.png'), kfold_rmseplot_bymodel_2wayfacet, height = 6, width = 7, dpi = 400)

ggsave('~/google_drive/NASABiodiversityWG/FlavorsOfDiversityPaper/figures_pdf/fig6.pdf', kfold_rmseplot_bymodel_2wayfacet, height = 6, width = 7)

# Table of fit stats ------------------------------------------------------

# By taxon, diversity level, flavor, model.
# Include RMSE, kfold RMSE, and R2 (with CIs)
# Also do a version with relative RMSEs to see if it's any better

fit_table <- all_rmse %>%
  mutate(taxon = factor(taxon, labels = c('birds','trees')),
         RMSE = paste0(round(RMSE_mean,2), ' [', round(RMSE_q025,2), ',', round(RMSE_q975,2), ']'),
         kfold_RMSE = paste0(round(kfold_RMSE_mean,2), ' [', round(kfold_RMSE_q025,2), ',', round(kfold_RMSE_q975,2), ']'),
         Rsquared = paste0(round(r2,2), ' [', round(r2_q025,2), ',', round(r2_q975,2), ']')) %>%
  select(taxon, rv, flavor, model, RMSE, kfold_RMSE, Rsquared) %>%
  arrange(taxon, rv, flavor, model)

write.csv(fit_table, file = '~/google_drive/NASABiodiversityWG/FlavorsOfDiversityPaper/supptable_fitstats.csv', row.names = FALSE)


fit_table_relative <- all_rmse %>%
  mutate(taxon = factor(taxon, labels = c('birds','trees')),
         RMSE = paste0(round(RMSE_mean_relative,2), ' [', round(RMSE_q025_relative,2), ',', round(RMSE_q975_relative,2), ']'),
         kfold_RMSE = paste0(round(kfold_RMSE_mean_relative,2), ' [', round(kfold_RMSE_q025_relative,2), ',', round(kfold_RMSE_q975_relative,2), ']'),
         Rsquared = paste0(round(r2,2), ' [', round(r2_q025,2), ',', round(r2_q975,2), ']')) %>%
  select(taxon, rv, flavor, model, RMSE, kfold_RMSE, Rsquared) %>%
  arrange(taxon, rv, flavor, model)

write.csv(fit_table_relative, file = '~/google_drive/NASABiodiversityWG/FlavorsOfDiversityPaper/supptable_fitstats_withrelativeRMSEs.csv', row.names = FALSE)

# WAIC table added 08 May 2019
waic_table <- model_waic %>%
  mutate(model = factor(model, levels=c('space','climate','geo','full'), labels=c('space only', 'space+climate','space+geodiversity','space+climate+geodiversity')),
         taxon = factor(taxon, labels = c('birds', 'trees')),
         WAIC = round(WAIC, 1),
         WAIC_SE = round(WAIC_SE, 1)) %>%
  select(taxon, rv, model, WAIC, WAIC_SE) %>%
  arrange(taxon, rv, model)

write.csv(waic_table, file = '~/google_drive/NASABiodiversityWG/FlavorsOfDiversityPaper/supptable_WAICs.csv', row.names = FALSE)
