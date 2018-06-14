# Make data needed to fit models for "flavors" analysis
# QDR NASABioXGeo 14 June 2018

# This new version, on 14 June, uses TRI (topographic ruggedness, mean neighbor pixel difference), instead of SD, for geodiversity.

# Load biodiversity and geodiversity data ---------------------------------

library(dplyr)
fp <- '/mnt/research/nasabio/data'

# Load adjacency matrices
load(file.path(fp, 'ecoregions/ecoregion_adjacency.RData'))

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
fiabio <- read.csv(file.path(fp, 'fia/biodiversity_CSVs/fia_bio_formixedmodels_50k.csv'), stringsAsFactors = FALSE) %>%
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
set.seed(38182)
fiasub_radius <- SRS_iterative_N1(fia_aea_noedge, radius = 20e3, n = 3000, point = sample(length(fia_aea_noedge),1), show_progress = TRUE)
fiaspsub <- fiasp_noedge[fiasub_radius, ]
fiabio <- subset(fiabio, PLT_CN %in% fiaspsub$CN)
fiageo <- subset(fiageo, PLT_CN %in% fiaspsub$CN)

# Unfortunately, get rid of the few plots where beta is 0 or 1
fiabio <- fiabio %>%
  mutate_at(vars(contains('beta_td')), function(x) if_else(x > 0 & x < 1, x, as.numeric(NA)))

# Save data so that models can be fit in parallel -------------------------

# Filter out the ones that are not in any ecoregion
# Edit 30 May: only do this for TNC
in_eco <- function(dat) dat$BCR %in% dimnames(bcr_bin)[[1]] & dat$TNC %in% dimnames(tnc_bin)[[1]] & dat$HUC4 %in% dimnames(huc_bin)[[1]]
in_eco <- function(dat) dat$TNC %in% dimnames(tnc_bin)[[1]]

fia_in_eco <- in_eco(fiageo)
bbs_in_eco <- in_eco(bbsgeo)

fiageo <- fiageo[fia_in_eco,]
fiabio <- fiabio[fia_in_eco,]
bbsgeo <- bbsgeo[bbs_in_eco,]
bbsbio <- bbsbio[bbs_in_eco,]

# 50 k version
save(fiageo, fiabio, bcr_bin, huc_bin, tnc_bin, file = '/mnt/research/nasabio/temp/fia_spatial_mm_dat_50k.RData')
save(bbsgeo, bbsbio, bcr_bin, huc_bin, tnc_bin, file = '/mnt/research/nasabio/temp/bbs_spatial_mm_dat_50k.RData')
