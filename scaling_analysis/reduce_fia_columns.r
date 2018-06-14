# Reduce FIA data to a reasonable number of columns for the scaling analysis
# QDR NASABioXGeo 14 June 2018

# Use native resolution, and use TRI (ruggedness)

######################## Geo

library(dplyr)
fp <- '/mnt/research/nasabio/data/fia'
fpgeo <- file.path(fp, 'geodiversity_CSVs')
bio1kmean <- read.csv(file.path(fpgeo, 'fia_bio1k_mean_wide.csv'), stringsAsFactors = FALSE)
elevmean <- read.csv(file.path(fpgeo, 'fia_elev_mean_wide.csv'), stringsAsFactors = FALSE)
othermean <- read.csv(file.path(fpgeo, 'fia_other_mean_wide.csv'), stringsAsFactors = FALSE)
othersd <- read.csv(file.path(fpgeo, 'fia_other_sd_wide.csv'), stringsAsFactors = FALSE)

fiageo <- Reduce(left_join,
				 list(
				 bio1kmean %>% 
					select(PLT_CN, contains('_5_'), contains('_20_'), contains('_100_')) %>%
					select(PLT_CN, contains('bio1_'), contains('bio12_')) %>%
					select(-contains('roughness')),
				 elevmean %>%
					select(PLT_CN, contains('_5_'), contains('_20_'), contains('_100_')) %>%
					select(PLT_CN, contains('elevation_30m')) %>%
					select(-contains('roughness')),
				 othermean %>%
					select(PLT_CN, contains('_5_'), contains('_20_'), contains('_100_')) %>%
					select(PLT_CN, contains('dhi_gpp'), contains('human_footprint')) %>%
					select(-contains('5k'), -contains('roughness')),
				 othersd %>% 
					select(PLT_CN, contains('_5_'), contains('_20_'), contains('_100_')) %>%
					select(PLT_CN, contains('geological_age'), contains('soil_type')) %>%
					select(-contains('diversity'), -contains('geological_age_5k')) ))
					
fiaeco <- read.csv(file.path(fp, 'fia_ecoregions.csv'))
fiageo <- right_join(fiaeco, fiageo)
write.csv(fiageo, file.path(fpgeo, 'fia_geo_forscalinganalysis.csv'), row.names = FALSE)	

######################### Bio
fpbio <- file.path(fp, 'biodiversity_CSVs')
ad <- read.csv(file.path(fpbio, 'fiausa_natural_alpha.csv'), stringsAsFactors = FALSE)
bd <- read.csv(file.path(fpbio, 'fiausa_natural_beta.csv'), stringsAsFactors = FALSE)
gd <- read.csv(file.path(fpbio, 'fiausa_natural_gamma.csv'), stringsAsFactors = FALSE)

rad <- 5
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
					
write.csv(fiabio, file.path(fpbio, 'fia_bio_forscalinganalysis.csv'), row.names = FALSE)