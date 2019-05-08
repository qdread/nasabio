# Create shaded polygon maps for mixed model fits
# This uses the new spatial mixed model coefficients. ***MULTIVARIATE!***
# QDR NASABIOXGEO 28 May 2018

# Edit 07 Jan 2019: update for new OS
# Edit 18 June: new predictor sets
# Edit 04 June: Also make maps of spatial effects only

task <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Define functions --------------------------------------------------------

# Make map (run fortify each time this function is called, slower but simpler)
# Return the plots as a list so that we can tile them in different ways
model_map <- function(coefs, fill_scale, regions, state_borders, bg_color = 'black', text_color = 'white', state_color = 'gray20', show_legend = TRUE) {
  # if significance is used, rename significance column so the code below doesn't break
  if ('significance' %in% names(coefs)) coefs <- rename(coefs, Estimate = significance)
  
  regions@data <- regions@data %>% 
    left_join(coefs)
  
  region_fort <- fortify(regions, region = 'id') %>% left_join(regions@data, by = 'id')
  
  blktheme <- theme_bw() + 
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          panel.grid = element_blank(), 
          panel.background = element_rect(color = bg_color, fill = bg_color), 
          panel.border = element_blank(), 
          plot.background = element_rect(fill = bg_color), 
          legend.position = c(0.15,0.1), 
          legend.key.width = unit(0.2, 'inches'),
          legend.direction = 'horizontal', 
          legend.title = element_blank(),
          legend.background = element_rect(fill = bg_color),
          legend.text = element_text(color = text_color))
  
  if (!show_legend) blktheme <- blktheme + theme(legend.position = 'none')
  
    ggplot(region_fort) +
      geom_polygon(aes(x=long, y=lat, group=group, fill=Estimate)) +
      geom_path(aes(x=long, y=lat, group=group), color = text_color, size = 0.25) +
      geom_path(data = state_borders, aes(x = long, y = lat, group = group), color = state_color) +
      fill_scale +
      coord_equal() +
      blktheme
}

arrangeMaps <- function(x, fpfig, prefix, titles, raw_names, text_color = 'white') {
  tw <- theme(plot.title = element_text(color = text_color))
  geo_name <- geo_names[which(prednames %in% x$parameter)] # For file name by geo variable
  x$bio_title <- titles[match(x$response, raw_names)] # Short title by bio variable
  x <- x[match(titles, x$bio_title),] # Put bio variables in correct order
  png(file.path(fpfig, paste0(prefix, '_', geo_name, '.png')), height = 9, width = 12, res = 400, units = 'in', type = 'cairo')
  grid.arrange(grobs = map2(x$maps, x$bio_title, function(p, name) ggplotGrob(p + ggtitle(name) + tw)), nrow = 3)
  dev.off()
  return('i just made a map :-)')
}

# Load region data --------------------------------------------------------

# Load coefficients
fpcoef <- '/mnt/research/nasabio/data/modelfits' # Cluster
fpregion <- '/mnt/research/nasabio/data/ecoregions'
fpstate <- '/mnt/home/qdr'


coef_all <- read.csv(file.path(fpcoef, 'multivariate_spatial_coef.csv'), stringsAsFactors = FALSE)

library(sp)
library(rgdal)
library(ggplot2)
library(dplyr)
library(purrr)
library(rgeos)
library(reshape2)

# Make a column for the parameter estimate and one that shows whether the CI is not zero.

coef_wide <- coef_all %>%
  filter(effect == 'coefficient', model == 'full') %>%
  select(-ecoregion, -effect, -rv, -model) %>%
  dcast(taxon + region + response + parameter ~ stat) %>%
  mutate(significance = case_when(
    Q2.5 > 0 ~ 'positive',
    Q97.5 < 0 ~ 'negative',
    TRUE ~ 'zero'
  ))

# Alternative version: spatial effects only (random, not including fixed effect)

spatialeff_wide <- coef_all %>%
  filter(effect == 'random', model == 'full') %>%
  select(-ecoregion, -effect, -rv, -model) %>%
  dcast(taxon + region + response + parameter ~ stat) %>%
  mutate(significance = case_when(
    Q2.5 > 0 ~ 'positive',
    Q97.5 < 0 ~ 'negative',
    TRUE ~ 'zero'
  ))



tnc <- readOGR(dsn = fpregion, layer = 'tnc_terr_ecoregions')
states <- read.csv(file.path(fpstate, 'states_albers.csv'), stringsAsFactors = FALSE)
load(file.path(fpstate, 'states_albers.RData'))

# Convert all to Albers with "region" in the name.
aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
tnc <- spTransform(tnc, CRSobj = CRS(aea_crs))

tnc@data <- tnc@data %>%
  mutate(id = rownames(tnc@data), region = as.character(ECODE_NAME))

# Subset out the regions that are outside the US.
tnc_unique <- unique(coef_all$region[coef_all$ecoregion == 'TNC'])
tnc <- subset(tnc, region %in% tnc_unique)

# Clip TNC to US boundaries
goodusabounds <- gUnaryUnion(states_albers)
tncdat <- tnc@data
tnc <- gIntersection(tnc, goodusabounds, byid = TRUE, id = row.names(tnc@data))
tnc <- SpatialPolygonsDataFrame(tnc, tncdat)

# Create all maps ------------------------------------------

rbfill <- scale_fill_gradient2(low = "#4575B4", high = "#D73027", midpoint = 0, limits = c(-3.6, 3.6), labels = -3:3, breaks = -3:3)
signif_fill <- scale_fill_manual(values = c(negative = "#4575B4", positive = "#D73027", zero = 'gray50'))

library(gridExtra)
fpfig <- '/mnt/research/nasabio/figs/mv_coefficient_maps'

bio_titles <- c('alpha taxonomic', 'alpha phylogenetic', 'alpha functional', 'beta taxonomic', 'beta phylogenetic', 'beta functional', 'gamma taxonomic', 'gamma phylogenetic', 'gamma functional')
bio_names <- c("alpha_richness", "alpha_phy_pa", "alpha_func_pa",
               "beta_td_sorensen_pa", "beta_phy_pa", "beta_func_pa",
               "gamma_richness", "gamma_phy_pa", "gamma_func_pa")
geo_names <- c('elevation_diversity','temperature_mean','geol_age_diversity','soil_diversity','precip_mean','gpp_sd', 'Intercept')
prednames <- c('elevation_5k_tri_50_mean', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'dhi_gpp_5k_tri_50_mean')

cols <- c(bg = 'white', text = 'black', state = 'gray20')

# Create list of maps for each predictor by response combo (9 predictors x 8 responses)

if (task == 1) {
  maps_bbs_coef <- coef_wide %>%
    filter(taxon == 'bbs') %>%
    dplyr::select(response, parameter, region, Estimate) %>%
    rename(ECODE_NAME = region) %>%
    group_by(response, parameter) %>%
    do(maps = model_map(., rbfill, tnc, states, bg_color=cols['bg'], text_color=cols['text'], state_color=cols['state']))
  maps_bbs_coef %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfig, prefix = 'bbs_coef', titles = bio_titles, raw_names = bio_names, text_color = cols['text']))
 }

if (task == 2) {
  maps_fia_coef <- coef_wide %>%
    filter(taxon == 'fia') %>%
    dplyr::select(response, parameter, region, Estimate) %>%
    rename(ECODE_NAME = region) %>%
    group_by(response, parameter) %>%
    do(maps = model_map(., rbfill, tnc, states, bg_color=cols['bg'], text_color=cols['text'], state_color=cols['state']))
  maps_fia_coef %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfig, prefix = 'fia_coef', titles = bio_titles, raw_names = bio_names, text_color = cols['text']))
}

if (task == 3) {
  maps_bbs_sig <- coef_wide %>%
    filter(taxon == 'bbs') %>%
    dplyr::select(response, parameter, region, significance) %>%
    rename(ECODE_NAME = region) %>%
    group_by(response, parameter) %>%
    do(maps = model_map(., signif_fill, tnc, states, bg_color=cols['bg'], text_color=cols['text'], state_color=cols['state'], show_legend = FALSE))
  maps_bbs_sig %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfig, prefix = 'bbs_signif', titles = bio_titles, raw_names = bio_names, text_color = cols['text']))
}

if (task == 4) {
  maps_fia_sig <- coef_wide %>%
    filter(taxon == 'fia') %>%
    dplyr::select(response, parameter, region, significance) %>%
    rename(ECODE_NAME = region) %>%
    group_by(response, parameter) %>%
    do(maps = model_map(., signif_fill, tnc, states, bg_color=cols['bg'], text_color=cols['text'], state_color=cols['state'], show_legend = FALSE))
  maps_fia_sig %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfig, prefix = 'fia_signif', titles = bio_titles, raw_names = bio_names, text_color = cols['text']))
}

if (task == 5) {
  maps_bbs_spat <- spatialeff_wide %>%
    filter(taxon == 'bbs') %>%
    dplyr::select(response, parameter, region, Estimate) %>%
    rename(ECODE_NAME = region) %>%
    group_by(response, parameter) %>%
    do(maps = model_map(., rbfill, tnc, states, bg_color=cols['bg'], text_color=cols['text'], state_color=cols['state']))
  maps_bbs_spat %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfig, prefix = 'bbs_spatial', titles = bio_titles, raw_names = bio_names, text_color = cols['text']))
}

if (task == 6) {
  maps_fia_spat <- spatialeff_wide %>%
    filter(taxon == 'fia') %>%
    dplyr::select(response, parameter, region, Estimate) %>%
    rename(ECODE_NAME = region) %>%
    group_by(response, parameter) %>%
    do(maps = model_map(., rbfill, tnc, states, bg_color=cols['bg'], text_color=cols['text'], state_color=cols['state']))
  maps_fia_spat %>%
    group_by(parameter) %>%
    do(foo = arrangeMaps(., fpfig = fpfig, prefix = 'fia_spatial', titles = bio_titles, raw_names = bio_names, text_color = cols['text']))
}
