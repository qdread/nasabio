# Create data frames of geodiversity and biodiversity data for the scaling analysis
# QDR Nasabioxgeo 14 June 2018

# Note that BBS must be the within-route biodiversity!

# This follows more or less the same procedure as multivariate_makedata.r

# Load biodiversity and geodiversity data ---------------------------------

library(dplyr)
fp <- '/mnt/research/nasabio/data'

# Load adjacency matrices
load(file.path(fp, 'ecoregions/ecoregion_adjacency.RData'))

# BBS
bbsbio <- read.csv(file.path(fp, 'bbs/biodiversity_CSVs/bbs_allbio_withinroute_wide.csv'), stringsAsFactors = FALSE)
bbsgeo <- read.csv(file.path(fp, 'bbs/geodiversity_CSVs/bbs_allgeo_wide.csv'), stringsAsFactors = FALSE)

# Use only 5 km radius for response variable (biodiv).
# For now keep all 9 of the variables, but we are probably only going to use alpha richness.
bbsbio <- bbsbio %>% 
  select(rteNo,
         alpha_richness_5, alpha_MPD_pa_z_5, alpha_MPDfunc_pa_z_5,
         beta_td_sorensen_pa_5, beta_pd_pairwise_pa_z_5, beta_fd_pairwise_pa_z_5,
         gamma_richness_5, gamma_MPD_pa_z_5, gamma_MPDfunc_pa_z_5) %>%
  rename(alpha_richness = alpha_richness_5, alpha_phy_pa = alpha_MPD_pa_z_5, alpha_func_pa = alpha_MPDfunc_pa_z_5,
         beta_td_sorensen_pa = beta_td_sorensen_pa_5, beta_phy_pa = beta_pd_pairwise_pa_z_5, beta_func_pa = beta_fd_pairwise_pa_z_5,
         gamma_richness = gamma_richness_5, gamma_phy_pa = gamma_MPD_pa_z_5, gamma_func_pa = gamma_MPDfunc_pa_z_5)
  

bbssp <- bbsgeo[,c('rteNo','lat','lon')]
bbsgeo <- bbsgeo %>%
  select(rteNo, HUC4, BCR, TNC, contains('_5_'), contains('_20_'), contains('_100_')) %>%
  select(-contains('roughness'), -contains('richness'), -contains('night'), -contains('mode'), -contains('sd')) %>%
  select(rteNo, HUC4, BCR, TNC,
         contains('bio1_1k'), contains('bio12_1k'), 
         contains('elevation_30m'), contains('soil'), contains('geological_age_1k'),
         contains('gpp_1k'), contains('footprint_1k')) %>%
  mutate(HUC4 = if_else(nchar(HUC4) == 3, paste0('0',HUC4), as.character(HUC4)))

# FIA (only needed variables are in these dfs, rest are on server)
# Use effective species number for Shannon
fiabio <- read.csv(file.path(fp, 'fia/biodiversity_CSVs/fia_bio_forscalinganalysis.csv'), stringsAsFactors = FALSE) %>%
  mutate(alpha_shannon = exp(alpha_shannon), gamma_shannon = exp(gamma_shannon)) %>%
  rename(alpha_effspn = alpha_shannon, gamma_effspn = gamma_shannon)
fiageo <- read.csv(file.path(fp, 'fia/geodiversity_CSVs/fia_geo_forscalinganalysis.csv'), stringsAsFactors = FALSE)
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

# Use 50km to Can/Mex. 
bbscoast <- flag_coast_plots(bbs_aea, radius = 50e3, border_countries = c('Canada','Mexico'))
noedge_rte <- bbssp$rteNo[!bbscoast$is_edge]
bbsbio <- subset(bbsbio, rteNo %in% noedge_rte)
bbsgeo <- subset(bbsgeo, rteNo %in% noedge_rte)

# Unfortunately, get rid of the few plots where beta is 0 or 1 (only 1 plot where it is 1, actually)
bbsbio <- bbsbio %>%
  mutate_at(vars(contains('beta_td')), function(x) if_else(x > 0 & x < 1, x, as.numeric(NA)))


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

# Filter out the ones that are not in a TNC ecoregion
in_eco <- function(dat) dat$TNC %in% dimnames(tnc_bin)[[1]]

fia_in_eco <- in_eco(fiageo)
bbs_in_eco <- in_eco(bbsgeo)

fiageo <- fiageo[fia_in_eco,]
fiabio <- fiabio[fia_in_eco,]
bbsgeo <- bbsgeo[bbs_in_eco,]
bbsbio <- bbsbio[bbs_in_eco,]

save(fiageo, fiabio, bcr_bin, huc_bin, tnc_bin, file = '/mnt/research/nasabio/temp/fia_scaling_dat.RData')
save(bbsgeo, bbsbio, bcr_bin, huc_bin, tnc_bin, file = '/mnt/research/nasabio/temp/bbs_scaling_dat.RData')
