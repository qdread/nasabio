# Create shaded polygon maps for mixed model fits
# This uses the new spatial mixed model coefficients.
# QDR NASABIOXGEO 27 Apr 2018


# Version created 30 Apr: Do in parallel on cluster
# Modified 1 May: Clip map to USA borders
# Modified 8 May: Do 50 km and 100 km separately, both with new data.
# Modified 10 May: plot the pop level coefficients, not the random effects
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
radius <- 50 # Or 100

# Define functions --------------------------------------------------------

# Make map (run fortify each time this function is called, slower but simpler)
# Return the plots as a list so that we can tile them in different ways
model_map <- function(coefs, fill_scale, regions, state_borders) {
  regions@data <- regions@data %>% 
    left_join(coefs)
  
  region_fort <- fortify(regions, region = 'id') %>% left_join(regions@data, by = 'id')
  
  blktheme <- theme_bw() + 
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          panel.grid = element_blank(), 
          panel.background = element_rect(color = 'black', fill = 'black'), 
          panel.border = element_blank(), 
          plot.background = element_rect(fill = 'black'), 
          legend.position = c(0.13,0.1), 
          legend.direction = 'horizontal', 
          legend.title = element_blank(),
          legend.background = element_rect(fill = 'black'),
          legend.text = element_text(color = 'white'))
  
    ggplot(region_fort) +
      geom_polygon(aes(x=long, y=lat, group=group, fill=Estimate)) +
      geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.25) +
      geom_path(data = state_borders, aes(x = long, y = lat, group = group), color = 'gray20') +
      fill_scale +
      coord_equal() +
      blktheme
}

arrangeMaps <- function(x, fpfig, region_name, titles, raw_names, div_type = 'incidence', rad = radius) {
  geo_name <- geo_names[which(prednames %in% x$parameter)] # For file name by geo variable
  x$bio_title <- titles[match(x$rv, raw_names)] # Short title by bio variable
  x <- x[match(titles, x$bio_title),] # Put bio variables in correct order
  png(file.path(fpfig, paste0(region_name, '_', rad, 'k_', div_type, '_', geo_name, '.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(x$maps, x$bio_title, function(p, name) ggplotGrob(p + ggtitle(name) + tw)), nrow = 3)
  dev.off()
  return('i just made a map :-)')
}

# Load region data --------------------------------------------------------

# Load coefficients
fpcoef <- '/mnt/research/nasabio/data/modelfits' # Cluster
fpregion <- '/mnt/research/nasabio/data/ecoregions'
fphuc <- fpregion
fpstate <- '~'


coef_all <- read.csv(file.path(fpcoef, paste0('spatial_coef_', radius, 'k.csv')), stringsAsFactors = FALSE)

# Edit HUC names
hucidx <- which(coef_all$ecoregion == 'HUC4' & nchar(coef_all$region) == 3)

coef_all$region[hucidx] <- paste0('0', coef_all$region[hucidx])

library(sp)
library(rgdal)
library(ggplot2)
library(dplyr)
library(purrr)
library(rgeos)

huc4 <- readOGR(dsn = fphuc, layer = 'HU4_CONUS_Alb')
bcr <- readOGR(dsn = fpregion, layer = 'BCR_Terrestrial_master')
tnc <- readOGR(dsn = fpregion, layer = 'tnc_terr_ecoregions')
states <- read.csv(file.path(fpstate, 'states_albers.csv'), stringsAsFactors = FALSE)
load(file.path(fpstate, 'states_albers.RData'))

# Convert all to Albers with "region" in the name.
aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
bcr <- spTransform(bcr, CRSobj = CRS(aea_crs))
tnc <- spTransform(tnc, CRSobj = CRS(aea_crs))

huc4@data <- huc4@data %>%
  mutate(id = rownames(huc4@data), region = as.character(HUC4))

bcr@data <- bcr@data %>%
  mutate(id = rownames(bcr@data), region = as.character(BCRNAME))

tnc@data <- tnc@data %>%
  mutate(id = rownames(tnc@data), region = as.character(ECODE_NAME))

# Subset out the regions that are outside the US.
bcr_unique <- unique(coef_all$region[coef_all$ecoregion == 'BCR'])
tnc_unique <- unique(coef_all$region[coef_all$ecoregion == 'TNC'])
bcr <- subset(bcr, region %in% bcr_unique & COUNTRY %in% 'USA' & !PROVINCE_S %in% 'ALASKA')
tnc <- subset(tnc, region %in% tnc_unique)

# Clip TNC to US boundaries
goodusabounds <- gUnaryUnion(states_albers)
tncdat <- tnc@data
tnc <- gIntersection(tnc, goodusabounds, byid = TRUE, id = row.names(tnc@data))
tnc <- SpatialPolygonsDataFrame(tnc, tncdat)

# Create all maps ------------------------------------------

rbfill <- scale_fill_gradient2(low = "#4575B4", high = "#D73027", midpoint = 0)

library(gridExtra)
fpbbs <- '/mnt/research/nasabio/figs/bbs_coefficient_maps'
if (radius == 100) prednames <- c('elevation_5k_100_sd', 'bio1_5k_100_mean', 'geological_age_5k_100_diversity', 'soil_type_5k_100_diversity', 'bio12_5k_100_mean', 'bio12_5k_100_sd', 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean')
if (radius == 50) prednames <- c('elevation_5k_50_sd', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'bio12_5k_50_sd', 'dhi_gpp_5k_50_sd', 'human_footprint_5k_50_mean')

bio_titles <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'beta PD', 'gamma PD', 'alpha FD', 'beta FD', 'gamma FD')
bio_names <- c("alpha_richness", "beta_td_sorensen_pa", "gamma_richness",
               "alpha_phy_pa", "beta_phy_pa", "gamma_phy_pa", 
               "alpha_func_pa", "beta_func_pa", "gamma_func_pa")
geo_names <- c('elevation_sd','temperature_mean','geol_age_diversity','soil_diversity','precip_mean','precip_sd','gpp_sd','footprint_mean')
tw <- theme(plot.title = element_text(color = 'white'))

# Separate incidence and abundance
fpfia <- '/mnt/research/nasabio/figs/fia_coefficient_maps'

# Index the incidence based and abundance based variables
# Note that there is no beta pd and fd for FIA
fia_bio_names_incid <- c("alpha_richness", "beta_td_sorensen_pa", "gamma_richness", 
                         "alpha_phy_pa", "beta_phy_pa", "gamma_phy_pa",
                         "alpha_func_pa", "beta_func_pa", "gamma_func_pa")
fia_bio_names_abund <- c("alpha_effspn", "beta_td_sorensen", "gamma_effspn",
                         "alpha_phy", "beta_phy", "gamma_phy",
                         "alpha_func", "beta_func", "gamma_func")

bio_titles_incidence <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'beta PD', 'gamma PD', 'alpha FD', 'beta FD', 'gamma FD')
bio_titles_abundance <- paste(bio_titles_incidence, 'abundance')


# Create list of maps for each predictor by response combo (9 predictors x 8 responses)
if (task == 1) {
  maps_bbs_huc <- coef_all %>%
    filter(taxon == 'bbs', effect == 'coefficient', ecoregion == 'HUC4') %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(HUC4 = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, huc4, states))
  maps_bbs_huc %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpbbs, region_name = 'HUC4', titles = bio_titles, raw_names = bio_names))
  
}

if (task == 2) {
  maps_bbs_bcr <- coef_all %>%
    filter(taxon == 'bbs', effect == 'coefficient', ecoregion == 'BCR') %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(BCRNAME = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, bcr, states))  
  maps_bbs_bcr %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpbbs, region_name = 'BCR', titles = bio_titles, raw_names = bio_names))
  
  
}

if (task == 3) {
  maps_bbs_tnc <- coef_all %>%
    filter(taxon == 'bbs', effect == 'coefficient', ecoregion == 'TNC') %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(ECODE_NAME = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, tnc, states))
  maps_bbs_tnc %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpbbs, region_name = 'TNC', titles = bio_titles, raw_names = bio_names))
  
  
}

if (task == 4) {
  maps_fia_huc_incid <- coef_all %>%
    filter(taxon == 'fia', effect == 'coefficient', ecoregion == 'HUC4', rv %in% fia_bio_names_incid) %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(HUC4 = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, huc4, states))
  maps_fia_huc_incid %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfia, region_name = 'HUC4', titles = bio_titles_incidence, raw_names = fia_bio_names_incid))
  
}

if (task == 5) {
  maps_fia_bcr_incid <- coef_all %>%
    filter(taxon == 'fia', effect == 'coefficient', ecoregion == 'BCR', rv %in% fia_bio_names_incid) %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(BCRNAME = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, bcr, states))  
  maps_fia_bcr_incid %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfia, region_name = 'BCR', titles = bio_titles_incidence, raw_names = fia_bio_names_incid))
}

if (task == 6) {
  maps_fia_tnc_incid <- coef_all %>%
    filter(taxon == 'fia', effect == 'coefficient', ecoregion == 'TNC', rv %in% fia_bio_names_incid) %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(ECODE_NAME = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, tnc, states))
  maps_fia_tnc_incid %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfia, region_name = 'TNC', titles = bio_titles_incidence, raw_names = fia_bio_names_incid))
}

if (task == 7) {
  maps_fia_huc_abund <- coef_all %>%
    filter(taxon == 'fia', effect == 'coefficient', ecoregion == 'HUC4', rv %in% fia_bio_names_abund) %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(HUC4 = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, huc4, states))
  maps_fia_huc_abund %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfia, region_name = 'HUC4', titles = bio_titles_abundance, raw_names = fia_bio_names_abund, div_type = 'abundance'))
}

if (task == 8) {
  maps_fia_bcr_abund <- coef_all %>%
    filter(taxon == 'fia', effect == 'coefficient', ecoregion == 'BCR', rv %in% fia_bio_names_abund) %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(BCRNAME = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, bcr, states))  
  maps_fia_bcr_abund %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfia, region_name = 'BCR', titles = bio_titles_abundance, raw_names = fia_bio_names_abund, div_type = 'abundance'))
}

if (task == 9) {
  maps_fia_tnc_abund <- coef_all %>%
    filter(taxon == 'fia', effect == 'coefficient', ecoregion == 'TNC', rv %in% fia_bio_names_abund) %>%
    dplyr::select(rv, parameter, region, Estimate) %>%
    rename(ECODE_NAME = region) %>%
    group_by(rv, parameter) %>%
    do(maps = model_map(., rbfill, tnc, states))
  maps_fia_tnc_abund %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfia, region_name = 'TNC', titles = bio_titles_abundance, raw_names = fia_bio_names_abund, div_type = 'abundance'))
}

