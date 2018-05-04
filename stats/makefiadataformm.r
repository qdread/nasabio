# Assemble data for mixed model for FIA (use same variables as BBS)
# Use only plots with plantations removed
# Use 100 km but later maybe use 50 km if this is too much.
# QDR 06 Mar 2018

# Edited 04 May 2018: Also make 50 km version.
# Edited 25 Apr 2018: Add phylogenetic and functional beta.

# Load biodiversity data

library(dplyr)
fp <- '/mnt/research/nasabio/data/fia'
ad <- read.csv(file.path(fp, 'fiausa_natural_alpha.csv'), stringsAsFactors = FALSE)
bd <- read.csv(file.path(fp, 'fiausa_natural_beta.csv'), stringsAsFactors = FALSE)
gd <- read.csv(file.path(fp, 'fiausa_natural_gamma.csv'), stringsAsFactors = FALSE)

rad <- 50

# Func, phy, and tax for alpha and gamma, but func and phy aren't ready yet for beta.
fiabio <- Reduce(left_join,
				 list(
				 ad %>% 
					filter(radius == rad) %>% 
					select(PLT_CN, richness, shannon, MPD_pa_z, MPD_z, MPDfunc_pa_z, MPDfunc_z) %>%
					rename(alpha_richness = richness, alpha_shannon = shannon, alpha_phy_pa = MPD_pa_z, alpha_phy = MPD_z, alpha_func_pa = MPDfunc_pa_z, alpha_func = MPDfunc_z),
				 bd %>%
					filter(radius == rad) %>%
					select(PLT_CN, beta_td_sorensen_pa, beta_td_sorensen, beta_pd_pairwise_pa_z, beta_pd_pairwise_z, beta_fd_pairwise_pa_z, beta_fd_pairwise_z) %>%
					rename(beta_phy_pa = beta_pd_pairwise_pa_z, beta_phy = beta_pd_pairwise_z, beta_func_pa = beta_fd_pairwise_pa_z, beta_func = beta_fd_pairwise_z),
				 gd %>% 
					filter(radius == rad) %>% 
					select(PLT_CN, richness, shannon, MPD_pa_z, MPD_z, MPDfunc_pa_z, MPDfunc_z) %>%
					rename(gamma_richness = richness, gamma_shannon = shannon, gamma_phy_pa = MPD_pa_z, gamma_phy = MPD_z, gamma_func_pa = MPDfunc_pa_z, gamma_func = MPDfunc_z)))
					
write.csv(fiabio, file.path(fp, 'fia_bio_formixedmodels_50k.csv'), row.names = FALSE)

# Load the geodiversity data frames that are needed: elevation, bio5k, geo and soil, dhi, footprint.

bio5kmean <- read.csv(file.path(fp, 'geodiv/fia_bio5k_mean_wide.csv'), stringsAsFactors = FALSE)
bio5ksd <- read.csv(file.path(fp, 'geodiv/fia_bio5k_sd_wide.csv'), stringsAsFactors = FALSE)
elevmean <- read.csv(file.path(fp, 'geodiv/fia_elev_mean_wide.csv'), stringsAsFactors = FALSE)
elevsd <- read.csv(file.path(fp, 'geodiv/fia_elev_sd_wide.csv'), stringsAsFactors = FALSE)
othermean <- read.csv(file.path(fp, 'geodiv/fia_other_mean_wide.csv'), stringsAsFactors = FALSE)
othersd <- read.csv(file.path(fp, 'geodiv/fia_other_sd_wide.csv'), stringsAsFactors = FALSE)

if (rad == 100) fiageo <- Reduce(left_join,
				 list(
				 bio5kmean %>% 
					select(PLT_CN, bio1_5k_100_mean, bio12_5k_100_mean),
				 bio5ksd %>% 
					select(PLT_CN, bio1_5k_100_sd, bio12_5k_100_sd),
				 elevmean %>%
					select(PLT_CN, elevation_5k_100_mean),
				 elevsd %>%
					select(PLT_CN, elevation_5k_100_sd),
				 othermean %>%
					select(PLT_CN, dhi_gpp_5k_100_mean, human_footprint_5k_100_mean),
				 othersd %>% 
					select(PLT_CN, dhi_gpp_5k_100_sd, human_footprint_5k_100_sd, geological_age_5k_100_diversity_geodiv, soil_type_5k_100_diversity_geodiv) ))
					
if (rad == 50) fiageo <- Reduce(left_join,
				 list(
				 bio5kmean %>% 
					select(PLT_CN, bio1_5k_50_mean, bio12_5k_50_mean),
				 bio5ksd %>% 
					select(PLT_CN, bio1_5k_50_sd, bio12_5k_50_sd),
				 elevmean %>%
					select(PLT_CN, elevation_5k_50_mean),
				 elevsd %>%
					select(PLT_CN, elevation_5k_50_sd),
				 othermean %>%
					select(PLT_CN, dhi_gpp_5k_50_mean, human_footprint_5k_50_mean),
				 othersd %>% 
					select(PLT_CN, dhi_gpp_5k_50_sd, human_footprint_5k_50_sd, geological_age_5k_50_diversity_geodiv, soil_type_5k_50_diversity_geodiv) ))
					
fiaeco <- read.csv(file.path(fp, 'fia_ecoregions.csv'))
fiageo <- right_join(fiaeco, fiageo)
write.csv(fiageo, file.path(fp, 'fia_geo_formixedmodels_50k.csv'), row.names = FALSE)	