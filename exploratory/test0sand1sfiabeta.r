# Script to locally run FIA beta-diversity model fits where zeroes and ones are simply removed, rather than adjusted as is done in the MS.


NC <- 2
NI <- 5000
NW <- 3000
delta <- 0.8

prednames <- c('elevation_5k_tri_50_mean', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'dhi_gpp_5k_tri_50_mean')
climate_prednames <- c('bio1_5k_50_mean', 'bio12_5k_50_mean')
geo_prednames <- c('elevation_5k_tri_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'dhi_gpp_5k_tri_50_mean')
alpha_resp <- c('alpha_richness', 'alpha_phy_pa', 'alpha_func_pa')
beta_resp <- c('beta_td_sorensen_pa', 'beta_phy_pa', 'beta_func_pa') 
gamma_resp <- c('gamma_richness', 'gamma_phy_pa', 'gamma_func_pa')

task_table <- expand.grid(taxon = c('fia','bbs'),
                          rv = c('alpha', 'beta', 'gamma'),
                          ecoregion = 'TNC',
                          model = c('full','climate','space', 'geo'),
                          fold = 0:63,
                          stringsAsFactors = FALSE)

taxon <- 'fia'
fold <- 0
rv <- beta_resp
# if(task_table$model[task] == 'climate') prednames <- climate_prednames
# if(task_table$model[task] == 'geo') prednames <- geo_prednames
# if(task_table$model[task] == 'space') prednames <- character(0)

ecoregion <- 'TNC'

source('~/Documents/GitHub/nasabio/stats/fit_mv_mm.r')

# Fit the model for the given response variable, taxon, and ecoregion
options(mc.cores = 3)

  load('~/Dropbox/projects/nasabiodiv/modelfits/fia_spatial_mm_dat_50k.RData')
  geodat <- fiageo
  biodat <- fiabio
  siteid <- 'PLT_CN'
  
  ### THIS IS THE MODIFIED PART FOR THIS EXPLORATORY ANALYSIS ###
  # Remove the (already converted to 0.001 and 0.999) values.
  biodat$beta_td_sorensen_pa[biodat$beta_td_sorensen_pa %in% range(biodat$beta_td_sorensen_pa, na.rm = TRUE)] <- NA
  
  # check how many are zeroes
  table(fiabio$beta_td_sorensen_pa == 0.001)
  table(fiabio$beta_td_sorensen_pa == 0.999)
  
  ### END MODIFICATION ###
  
  
  # Added 14 May: logit transform beta td.
  biodat$beta_td_sorensen_pa <- qlogis(biodat$beta_td_sorensen_pa)
  biodat$beta_td_sorensen <- qlogis(biodat$beta_td_sorensen)

  
  

# The following six ecoregions should not be used in any model fitting because they have too few data points. 
# They are primarily in Canada or Mexico with only a small portion of area in the USA, once buffer is deducted

exclude_regions <- c('NA0801', 'NA0808', 'NA0417', 'NA0514', 'NA1202', 'NA1301')

# Set data from the holdout set to missing, if task was designated as a k-fold task
# For "leave one region out" cross-validation, we just need to get rid of a single region for each fold

# Added 02 May 2019: include the ecoregion folds, less the excluded ones
fold_df <- read.csv('~/Dropbox/projects/nasabiodiv/ecoregion_folds.csv', stringsAsFactors = FALSE)
region_folds <- fold_df$TNC
region_folds <- region_folds[!grepl(paste(exclude_regions, collapse = '|'), region_folds)]

library(dplyr)

if (fold != 0) {
	# Join response variable data with the region ID, then set the appropriate values to NA
	biodat <- biodat %>% left_join(geodat[, c(siteid, 'TNC')])
	biodat$missing <- biodat$TNC == region_folds[fold]
}

# Modified 14 May: model all with Gaussian
distrib <- 'gaussian'
         
# Priors (added May 29)
# --------------------

# Edit 04 Jan 2019: temporarily remove all priors (add some back in on 05 Jan)
# Edit May 31: Add priors for FIA intercepts and for BBS alpha sdcar
# Edit June 14: Add sdcar priors and intercept priors on FIA beta, sd car priors on BBS beta
library(brms)
# 1st arg is df, 2nd is mu, 3rd is sigma for student t distribution
added_priors <- NULL
# if (task_table$rv[task] == 'alpha' & taxon == 'fia') {
#   added_priors <- c(set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alpharichness'),
# 					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alphaphypa'),
# 					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alphafuncpa') )
# } 
# if (task_table$rv[task] == 'beta' & taxon == 'fia') {
  added_priors <- c(set_prior('lognormal(1, 2)', class = 'sdcar', resp = 'betatdsorensenpa'),
					set_prior('lognormal(1, 2)', class = 'sdcar', resp = 'betaphypa'),
					set_prior('lognormal(1, 2)', class = 'sdcar', resp = 'betafuncpa'),
					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betatdsorensenpa'),
					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betaphypa'),
					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betafuncpa')					)
# }
# if (task_table$rv[task] == 'beta' & taxon == 'fia') {
#   added_priors <- c(set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betatdsorensenpa'),
# 					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betaphypa'),
# 					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betafuncpa') )					
# }
# if (task_table$rv[task] == 'gamma' & taxon == 'fia') {
#   added_priors <- c(set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'gammarichness'),
# 					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'gammaphypa'),
# 					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'gammafuncpa') )
# } 
# if (task_table$rv[task] == 'alpha' & taxon == 'bbs') {
#   added_priors <- c(set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'alpharichness'),
# 					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'alphaphypa'),
# 					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'alphafuncpa'),
# 					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alpharichness'),
# 					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alphaphypa'),
# 					set_prior('student_t(5, 0, 2)', class = 'Intercept', resp = 'alphafuncpa') )
# }
# if (task_table$rv[task] == 'beta' & taxon == 'bbs') {
#   added_priors <- c(set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betatdsorensenpa'),
# 					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betaphypa'),
# 					set_prior('lognormal(1, 1)', class = 'sdcar', resp = 'betafuncpa'),
# 					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betatdsorensenpa'),
# 					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betaphypa'),
# 					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'betafuncpa') )
# }
# if (task_table$rv[task] == 'gamma' & taxon == 'bbs') {
#   added_priors <- c(set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'gammarichness'),
# 					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'gammaphypa'),
# 					set_prior('student_t(10, 0, 1)', class = 'Intercept', resp = 'gammafuncpa') )
# } 
# 			 
# --------------------				  
				  
if (ecoregion == 'HUC4') eco_mat <- huc_bin
if (ecoregion == 'BCR') eco_mat <- bcr_bin
if (ecoregion == 'TNC') eco_mat <- tnc_bin

fit <- fit_mv_mm(pred_df = geodat, 
                 resp_df = biodat, 
                 pred_vars = prednames, 
                 resp_vars = rv, 
                 id_var = siteid, 
                 region_var = ecoregion, 
                 distribution = distrib, 
                 adj_matrix = eco_mat,
                 priors = added_priors,
                 n_chains = 3,
                 n_iter = 5000,
                 n_warmup = 3000,
                 delta = delta,
                 missing_data = fold > 0,
                 exclude_locations = exclude_regions
)

# Save all fits -- do in scratch space because so big
# save(fit, file = paste0('/mnt/gs18/scratch/groups/nasabio/modelfits/fit',task,'.RData'))

# Save the fit temporarily.
save(fit, file = '~/Dropbox/projects/nasabiodiv/modelfits/testfiafitremove0and1.RData')

# Extract fixed effect coefficients from model fit.

library(tidyr)
library(ggplot2)

model_coef <- fit$coef

model_fixef <- model_coef %>% 
  filter(effect %in% 'fixed') %>%
  separate(parameter, into = c('response', 'parameter'), sep= '_', extra = 'merge') %>%
  select(response, parameter, stat, value) %>%
  group_by(response, parameter) %>%
  spread(stat, value)

ggplot(model_fixef %>% filter(!parameter %in% 'Intercept'), aes(x = parameter, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  facet_wrap(~ response) +
  geom_pointrange() +
  coord_flip()

# Load the non additive one
model_coef_orig <- read.csv('~/Dropbox/projects/nasabiodiv/modelfits/multivariate_spatial_coef.csv', stringsAsFactors = FALSE)

model_fixef_orig <- model_coef_orig %>%
  filter(effect %in% 'fixed', model %in% 'full', taxon %in% 'bbs', rv %in% 'beta') %>%
  mutate(response = gsub('_', '', response)) %>%
  select(response, parameter, stat, value) %>%
  group_by(response, parameter) %>%
  spread(stat, value)

model_fixef_all <- rbind(data.frame(model = 'old', model_fixef_orig), data.frame(model = 'new', model_fixef)) %>%
  mutate(schnignificant = (Q2.5>0 & Q97.5>0) | (Q2.5<0 & Q97.5<0))

ggplot(model_fixef_all %>% filter(!parameter %in% 'Intercept'), aes(x = parameter, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  facet_grid(model ~ response) +
  geom_pointrange() +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'red') +
  coord_flip() +
  theme_bw()

# Make a slightly nicer figure so that it can be shown to the co-authors
param_labels <- c('climate: temp mean', 'climate: precip mean', 'geodiv: GPP', 'geodiv: elevation', 'geodiv: geological age', 'geodiv: soil type')
beta_labels <- c('betafuncpa' = 'Functional beta', 'betaphypa' = 'Phylogenetic beta', 'betatdsorensenpa' = 'Taxonomic beta')

ggplot(model_fixef_all %>% filter(!parameter %in% 'Intercept'), aes(x = parameter, y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  facet_wrap(~ response, labeller = labeller(response = beta_labels)) +
  geom_pointrange(aes(color = model)) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'blue') +
  scale_x_discrete(labels = param_labels) +
  coord_flip() +
  theme_bw() +
  scale_color_manual(values = c('blue', 'red'), labels = c('0 and 1 adjusted\n(used in MS)', '0 and 1 removed')) +
  theme(legend.position = 'bottom')
