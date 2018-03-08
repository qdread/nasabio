# Assemble data for mixed model for FIA (use same variables as BBS)
# Use only plots with plantations removed
# Use 100 km but later maybe use 50 km if this is too much.
# QDR 06 Mar 2018

# Load biodiversity data

library(dplyr)
fp <- '/mnt/research/nasabio/data/fia'
ad <- read.csv(file.path(fp, 'fiausa_natural_alpha.csv'), stringsAsFactors = FALSE)
bd <- read.csv(file.path(fp, 'fiausa_natural_betatd.csv'), stringsAsFactors = FALSE)
gd <- read.csv(file.path(fp, 'fiausa_natural_gamma.csv'), stringsAsFactors = FALSE)

# Func, phy, and tax for alpha and gamma, but func and phy aren't ready yet for beta.
fiabio <- Reduce(left_join,
				 list(
				 ad %>% 
					filter(radius == 100) %>% 
					select(PLT_CN, richness, shannon, MPD_pa_z, MPD_z, MPDfunc_pa_z, MPDfunc_z) %>%
					rename(alpha_richness = richness, alpha_shannon = shannon, alpha_phy_pa = MPD_pa_z, alpha_phy = MPD_z, alpha_func_pa = MPDfunc_pa_z, alpha_func = MPDfunc_z),
				 bd %>%
					filter(radius == 100) %>%
					select(PLT_CN, beta_td_sorensen_pa, beta_td_sorensen),
				 gd %>% 
					filter(radius == 100) %>% 
					select(PLT_CN, richness, shannon, MPD_pa_z, MPD_z, MPDfunc_pa_z, MPDfunc_z) %>%
					rename(gamma_richness = richness, gamma_shannon = shannon, gamma_phy_pa = MPD_pa_z, gamma_phy = MPD_z, gamma_func_pa = MPDfunc_pa_z, gamma_func = MPDfunc_z)))
					
write.csv(fiabio, file.path(fp, 'fia_bio_formixedmodels.csv'), row.names = FALSE)

# Load the geodiversity data frames that are needed: elevation, bio5k, geo and soil, dhi, footprint.

bio5kmean <- read.csv(file.path(fp, 'geodiv/fia_bio5k_mean_wide.csv'), stringsAsFactors = FALSE)
bio5ksd <- read.csv(file.path(fp, 'geodiv/fia_bio5k_sd_wide.csv'), stringsAsFactors = FALSE)
elevmean <- read.csv(file.path(fp, 'geodiv/fia_elev_mean_wide.csv'), stringsAsFactors = FALSE)
elevsd <- read.csv(file.path(fp, 'geodiv/fia_elev_sd_wide.csv'), stringsAsFactors = FALSE)
othermean <- read.csv(file.path(fp, 'geodiv/fia_other_mean_wide.csv'), stringsAsFactors = FALSE)
othersd <- read.csv(file.path(fp, 'geodiv/fia_other_sd_wide.csv'), stringsAsFactors = FALSE)

fiageo <- Reduce(left_join,
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
					
write.csv(fiageo, file.path(fp, 'fia_geo_formixedmodels.csv'), row.names = FALSE)

# Added 07 March: add ecoregions to this data frames
fiageo <- read.csv(file.path(fp, 'fia_geo_formixedmodels.csv'))
fiaeco <- read.csv(file.path(fp, 'fia_ecoregions.csv'))
fiageo <- right_join(fiaeco, fiageo)
write.csv(fiageo, file.path(fp, 'fia_geo_formixedmodels.csv'), row.names = FALSE)	