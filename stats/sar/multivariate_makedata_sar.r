# Make data needed to fit models for "flavors" analysis
# QDR NASABioXGeo 14 June 2018

# This new version uses TRI (topographic ruggedness, mean neighbor pixel difference), instead of SD, for geodiversity.
# Edited 04 Jan 2019: Change how 0 and 1 values are treated for FIA
# New version created 29 Apr 2019: Make a "real" distance matrix so that we can fit a SAR not a CAR.

# Load biodiversity and geodiversity data ---------------------------------

library(dplyr)
fp <- '/mnt/research/nasabio/data'

# BBS
bbsbio <- read.csv(file.path(fp, 'bbs/biodiversity_CSVs/bbs_allbio_wide.csv'), stringsAsFactors = FALSE)
bbsgeo <- read.csv(file.path(fp, 'bbs/geodiversity_CSVs/bbs_allgeo_wide.csv'), stringsAsFactors = FALSE)

# 50 km version -----------------------------------------------------------

# Use only 50 km radius.
bbsbio <- bbsbio %>% 
  select(rteNo,
         alpha_richness_50, alpha_MPD_pa_z_50, alpha_MPDfunc_pa_z_50,
         beta_td_sorensen_pa_50, beta_pd_pairwise_pa_z_50, beta_fd_pairwise_pa_z_50,
         gamma_richness_50, gamma_MPD_pa_z_50, gamma_MPDfunc_pa_z_50) %>%
  rename(alpha_richness = alpha_richness_50, alpha_phy_pa = alpha_MPD_pa_z_50, alpha_func_pa = alpha_MPDfunc_pa_z_50,
         beta_td_sorensen_pa = beta_td_sorensen_pa_50, beta_phy_pa = beta_pd_pairwise_pa_z_50, beta_func_pa = beta_fd_pairwise_pa_z_50,
         gamma_richness = gamma_richness_50, gamma_phy_pa = gamma_MPD_pa_z_50, gamma_func_pa = gamma_MPDfunc_pa_z_50)


bbssp <- bbsgeo[,c('rteNo','lat','lon')]
bbsgeo <- bbsgeo %>%
  select(rteNo, HUC4, BCR, TNC, contains('point'), matches('5k.*_50_')) %>%
  select(-contains('roughness'), -contains('richness'), -contains('night')) %>%
  select(rteNo, HUC4, BCR, TNC,
         contains('bio1_'), contains('bio12_'), contains('biocloud1_'),
         contains('elevation'), contains('soil'), contains('geological'),
         contains('footprint'), contains('gpp')) %>%
  select(-contains('point')) %>%
  mutate(HUC4 = if_else(nchar(HUC4) == 3, paste0('0',HUC4), as.character(HUC4)))

# FIA (only needed variables are in these dfs, rest are on server)
# Refer to makefiadataformm.r for making these reduced data frames from the huge data frames.
# Use effective species number for Shannon
fiabio <- read.csv(file.path(fp, 'fia/biodiversity_CSVs/updated_nov2018/fia_bio_formixedmodels_50k.csv'), stringsAsFactors = FALSE) %>%
  mutate(alpha_shannon = exp(alpha_shannon), gamma_shannon = exp(gamma_shannon)) %>%
  rename(alpha_effspn = alpha_shannon, gamma_effspn = gamma_shannon)
fiageo <- read.csv(file.path(fp, 'fia/geodiversity_CSVs/fia_geo_formixedmodels_50k.csv'), stringsAsFactors = FALSE)
fiageo <- fiageo %>%
  setNames(gsub('_geodiv','',names(.))) %>%
  mutate(HUC4 = if_else(nchar(HUC4) == 3, paste0('0',HUC4), as.character(HUC4)))


# Get rid of coasts and thin down plots -----------------------------------

library(sp)

aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
bbssp <- subset(bbssp, rteNo %in% bbsbio$rteNo)
bbs_aea <- SpatialPoints(coords=bbssp[,c('lon','lat')], proj4string = CRS('+proj=longlat +ellps=WGS84 +no_defs'))
bbs_aea <- spTransform(bbs_aea, CRS(aea_crs))

# Function for iterative search.
source('/mnt/research/nasabio/code/SRS_iterative.r')
# Functions for flagging edge plots
source('/mnt/research/nasabio/code/spatial_fns.r')

# Use 50km to Can/Mex. Don't thin BBS plots
bbscoast <- flag_coast_plots(bbs_aea, radius = 50e3, border_countries = c('Canada','Mexico'))
noedge_rte <- bbssp$rteNo[!bbscoast$is_edge]
bbsbio <- subset(bbsbio, rteNo %in% noedge_rte)
bbsgeo <- subset(bbsgeo, rteNo %in% noedge_rte)

# For FIA, use 10 km thinning and get rid of 50km to Can/Mex. This will reduce dataset to very manageable <10k plots.
fiasp <- read.csv('~/data/allfia.csv')
aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
fiasp <- subset(fiasp, CN %in% fiabio$PLT_CN)
fia_aea <- SpatialPoints(coords=fiasp[,c('ACTUAL_LON','ACTUAL_LAT')], proj4string = CRS('+proj=longlat +ellps=WGS84 +no_defs'))
fia_aea <- spTransform(fia_aea, CRS(aea_crs))

fiacoast <- flag_coast_plots(fia_aea, radius = 50e3, border_countries = c('Canada','Mexico'))

fia_aea_noedge <- fia_aea[!fiacoast$is_edge]
fiasp_noedge <- fiasp[!fiacoast$is_edge, ]
# Thin FIA plots to minimum 10 km separation
# Edited 30 May: 20 km instead, to try to equalize with BBS's sample size (ends up being ~ 3000 also)
set.seed(111)
fiasub_radius <- SRS_iterative_N1(fia_aea_noedge, radius = 20e3, n = 3000, point = sample(length(fia_aea_noedge),1), show_progress = TRUE)
fiaspsub <- fiasp_noedge[fiasub_radius, ]
fiabio <- subset(fiabio, PLT_CN %in% fiaspsub$CN)
fiageo <- subset(fiageo, PLT_CN %in% fiaspsub$CN)

# Added 4 Jan 2019: Instead of getting rid of the 0 and 1 plots, use a correction factor
fiabio <- fiabio %>%
  mutate_at(vars(contains('beta_td')), function(x) case_when(
	x == 0 ~ 0.001,
	x == 1 ~ 0.999,
	TRUE ~ x))

# Create distance matrix --------------------------------------------------

# This is now needed for the SAR.
# Make sure that the thinned out locations and the edge locations have been removed

# Subset projected object by the removed points
bbs_aea_noedge <- bbs_aea[!bbscoast$is_edge]
bbsdist <- spDists(bbs_aea_noedge)

fia_aea_noedge_thinned <- fia_aea_noedge[fiasub_radius]
fiadist <- spDists(fia_aea_noedge_thinned)

# Add custom k-fold column to geo data ------------------------------------

# (Don't do this for now. It can be added using the ecoregion groups.)
	
# Save data so that models can be fit in parallel -------------------------

# 50 k version
save(fiageo, fiabio, fiadist, file = '/mnt/research/nasabio/temp/fia_spatial_mm_dat_50k_sar.RData')
save(bbsgeo, bbsbio, bbsdist, file = '/mnt/research/nasabio/temp/bbs_spatial_mm_dat_50k_sar.RData')
