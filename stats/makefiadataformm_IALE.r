# Assemble data for mixed model for FIA (use same variables as BBS)
# Use only plots with plantations removed
# Use 100 km but later maybe use 50 km if this is too much.
# QDR 06 Mar 2018
# New version 05 Apr: IALE.

# Load biodiversity data

library(dplyr)
fp <- '/mnt/research/nasabio/data/fia'
ad <- read.csv(file.path(fp, 'fiausa_natural_alpha.csv'), stringsAsFactors = FALSE)
bd <- read.csv(file.path(fp, 'fiausa_natural_betatd.csv'), stringsAsFactors = FALSE)
gd <- read.csv(file.path(fp, 'fiausa_natural_gamma.csv'), stringsAsFactors = FALSE)

radii <- c(5, 10, 20, 50, 100, 200)

fiabio <- Reduce(left_join,
				 list(
				 ad %>% 
					filter(radius %in% radii) %>% 
					select(PLT_CN, radius, richness) %>%
					rename(alpha_richness = richness),
				 bd %>%
					filter(radius %in% radii) %>%
					select(PLT_CN, radius, beta_td_sorensen_pa),
				 gd %>% 
					filter(radius %in% radii) %>% 
					select(PLT_CN, radius, richness) %>%
					rename(gamma_richness = richness)))
					
write.csv(fiabio, '/mnt/research/nasabio/temp/fia_bio_formixedmodels_IALE.csv', row.names = FALSE)

bio5kmean <- read.csv(file.path(fp, 'geodiv/fia_bio5k_mean_wide.csv'), stringsAsFactors = FALSE)
bio5ksd <- read.csv(file.path(fp, 'geodiv/fia_bio5k_sd_wide.csv'), stringsAsFactors = FALSE)
elevmean <- read.csv(file.path(fp, 'geodiv/fia_elev_mean_wide.csv'), stringsAsFactors = FALSE)
elevsd <- read.csv(file.path(fp, 'geodiv/fia_elev_sd_wide.csv'), stringsAsFactors = FALSE)
othermean <- read.csv(file.path(fp, 'geodiv/fia_other_mean_wide.csv'), stringsAsFactors = FALSE)
othersd <- read.csv(file.path(fp, 'geodiv/fia_other_sd_wide.csv'), stringsAsFactors = FALSE)

# For IALE, need 6 radii, 5k resolution. From mean df, need roughness mean and tri mean, from sd df, need sd.
# Variables are bio12, elevation, footprint

radii <- paste0('_', c(5,10,20,50,100,200), '_')

fiageo <- Reduce(left_join,
				 list(
				 bio5kmean %>% 
					select(PLT_CN, contains('bio12_')) %>%
					select(PLT_CN, contains('_5k_')) %>%
					select(PLT_CN, contains('roughness'), contains('tri')) %>%
					select(PLT_CN, matches(paste0(radii, collapse = '|'))),
				 bio5ksd %>% 
					select(PLT_CN, contains('bio12_')) %>%
					select(PLT_CN, contains('_5k_')) %>%
					select(PLT_CN, contains('sd')) %>%
					select(PLT_CN, matches(paste0(radii, collapse = '|'))),
				 elevmean %>%
					select(PLT_CN, contains('elevation_')) %>%
					select(PLT_CN, contains('_5k_')) %>%
					select(PLT_CN, contains('roughness'), contains('tri')) %>%
					select(PLT_CN, matches(paste0(radii, collapse = '|'))),
				 elevsd %>%
					select(PLT_CN, contains('elevation_')) %>%
					select(PLT_CN, contains('_5k_')) %>%
					select(PLT_CN, contains('sd')) %>%
					select(PLT_CN, matches(paste0(radii, collapse = '|'))),
				 othermean %>%
					select(PLT_CN, contains('footprint')) %>%
					select(PLT_CN, contains('_5k_')) %>%
					select(PLT_CN, contains('roughness'), contains('tri')) %>%
					select(PLT_CN, matches(paste0(radii, collapse = '|'))),
				 othersd %>% 
					select(PLT_CN, contains('footprint')) %>%
					select(PLT_CN, contains('_5k_')) %>%
					select(PLT_CN, contains('sd')) %>%
					select(PLT_CN, matches(paste0(radii, collapse = '|'))) 
				))
					
fiaeco <- read.csv(file.path(fp, 'fia_ecoregions.csv'))
fiageo <- right_join(fiaeco[,c('PLT_CN','HUC4')], fiageo)
write.csv(fiageo, '/mnt/research/nasabio/temp/fia_geo_formixedmodels_IALE.csv', row.names = FALSE)
