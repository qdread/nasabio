# Create data frames for FIA and BBS variable selection
# Use the TRI variables, standard deviation variables, and mean variables

# Bind together data from BBS and FIA for variable selection
# ==========================================================

library(tidyverse)
fp <- '/mnt/research/nasabio/data'

bbsgeo <- read_csv(file.path(fp, 'bbs/geodiversity_CSVs/bbs_allgeo_wide.csv'))

fpgeo <- file.path(fp, 'fia/geodiversity_CSVs')
f <- function(x) x %>% select(PLT_CN, matches('5k.*_50_'))	

bio5kmean <- read_csv(file.path(fpgeo, 'fia_bio5k_mean_wide.csv')) %>% f
bio5ksd <- read_csv(file.path(fpgeo, 'fia_bio5k_sd_wide.csv')) %>% f
elevmean <- read_csv(file.path(fpgeo, 'fia_elev_mean_wide.csv')) %>% f
elevsd <- read_csv(file.path(fpgeo, 'fia_elev_sd_wide.csv')) %>% f
othermean <- read_csv(file.path(fpgeo, 'fia_other_mean_wide.csv')) %>% f
othersd <- read_csv(file.path(fpgeo, 'fia_other_sd_wide.csv')) %>% f

bbsgeo <- bbsgeo %>%
	select(rteNo, matches('5k.*_50_'))

fiageo <- Reduce(left_join, list(bio5kmean, bio5ksd, elevmean, elevsd, othermean, othersd)) %>% 
	setNames(gsub('_geodiv','',names(.)))
	
bbsgeo <- bbsgeo %>% rename(id = rteNo)
fiageo <- fiageo %>% rename(id = PLT_CN)

allgeo <- bbsgeo %>%
	select(-matches('biocloud')) %>%
	bind_rows(fiageo)
	
write.csv(allgeo, '/mnt/research/nasabio/temp/dat_for_variable_selection.csv', row.names = FALSE)


# Reload bound data and do some variable selection things
# =======================================================

library(tidyverse)
allgeo <- read_csv('/mnt/research/nasabio/temp/dat_for_variable_selection.csv')

# Don't use the roughness variables
allgeo_comp <- allgeo %>% 
	select(-contains('roughness')) %>%
	filter(complete.cases(.))

geopca <- prcomp(allgeo_comp[,-1], center = TRUE, scale = TRUE)

# What are the top few variables for each column?

for (i in 1:10) {
	topfive <- names(sort(abs(geopca$rotation[,i]), decreasing = TRUE))[1:5]
	print(geopca$rotation[topfive, i])
}

climatevars <- allgeo_comp %>%
	select(matches('bio.*_5k_50_mean'))
	
climatepca <- prcomp(climatevars, center = TRUE, scale = TRUE)
