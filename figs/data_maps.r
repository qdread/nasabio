# Maps to show raw data for NASABIOXGEO MS
# QDR/29 May 2018

# Edit 18 Jun 2018: new predictor variables.
# Edit 4 Jun 2018: change midpoint to gray, change labels

# Make map (run fortify each time this function is called, slower but simpler)
# Return the plots as a list so that we can tile them in different ways
point_map <- function(dat, color_scale, regions, state_borders, bg_color = 'black', text_color = 'white', state_color = 'gray20', fill_color = 'white', show_legend = TRUE) {
  
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
          legend.text = element_text(color = text_color, size = 6))
  
  if (!show_legend) blktheme <- blktheme + theme(legend.position = 'none')
  
  ggplot(region_fort) +
    geom_polygon(aes(x=long, y=lat, group=group), fill = fill_color) +
    geom_point(data = dat, aes(x = lon_aea, y = lat_aea, color = value), size = 1.25) +
    geom_path(aes(x=long, y=lat, group=group), color = text_color, size = 0.25) +
    geom_path(data = state_borders, aes(x = long, y = lat, group = group), color = state_color) +
    color_scale +
    coord_equal() +
    blktheme
}

library(dplyr)
library(reshape2)
library(ggplot2)
library(rgdal)
library(rgeos)
library(purrr)
library(gridExtra)

fpregion <- '/mnt/research/nasabio/data/ecoregions'
fpstate <- '~'
fpfig <- '/mnt/research/nasabio/figs/descriptivemaps'

rbcol <- scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(9,'RdYlBu')))

prednames <- c('elevation_5k_tri_50_mean', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'dhi_gpp_5k_tri_50_mean')

bio_titles <- c('alpha taxonomic', 'alpha phylogenetic', 'alpha functional', 'beta taxonomic', 'beta phylogenetic', 'beta functional', 'gamma taxonomic', 'gamma phylogenetic', 'gamma functional')
bio_names <- c("alpha_richness", "alpha_phy_pa", "alpha_func_pa",
               "beta_td_sorensen_pa", "beta_phy_pa", "beta_func_pa",
               "gamma_richness", "gamma_phy_pa", "gamma_func_pa")
geo_names <- c('elevation_diversity','temperature_mean','geol_age_diversity','soil_diversity','precip_mean','gpp_sd', 'intercept')

# Load data and coordinates
load('/mnt/research/nasabio/temp/bbs_spatial_mm_dat_50k.RData')
# Added April 30: Correction for outliers on beta functional
bbsbio$beta_func_pa[bbsbio$beta_func_pa < -10] <- NA
load('/mnt/research/nasabio/temp/fia_spatial_mm_dat_50k.RData')

# Load coordinates (fuzzed coordinates are OK for FIA)
# They are both in Albers
bbscoords <- read.csv('/mnt/research/nasabio/data/bbs/bbs_route_midpoints.csv')
fiacoords <- read.csv('/mnt/research/nasabio/data/fia/fia_fuzzed_coords.csv')

# Combine all the predictors (geo) into single data frame.
bbsgeo <- bbsgeo %>% left_join(bbscoords) %>% dplyr::select(-rteNo)
fiageo <- fiageo %>% left_join(fiacoords) %>% dplyr::select(-PLT_CN) %>%
  setNames(gsub('_fuzzed', '', names(.)))
allgeo <- bind_rows(bbsgeo, fiageo)

# Join coordinates with bio
bbsbio <- bbsbio %>% left_join(bbscoords) %>% filter(lat < 50)
fiabio <- fiabio %>% left_join(fiacoords) %>%
  setNames(gsub('_fuzzed', '', names(.)))

# Convert everything to longform
allgeo_long <- allgeo %>%
  select(-lon, -lat, -HUC4, -BCR, -TNC) %>%
  melt(id.vars = c('lon_aea', 'lat_aea')) %>%
  filter(variable %in% prednames)

bbsbio_long <- bbsbio %>%
  select(-lon, -lat) %>%
  melt(id.vars = c('rteNo', 'lon_aea', 'lat_aea')) %>%
  filter(variable %in% bio_names)

fiabio_long <- fiabio %>%
  select(-lon, -lat) %>%
  melt(id.vars = c('PLT_CN', 'lon_aea', 'lat_aea')) %>%
  filter(variable %in% bio_names)

# Load TNC info
tnc <- readOGR(dsn = fpregion, layer = 'tnc_terr_ecoregions')
states <- read.csv(file.path(fpstate, 'states_albers.csv'), stringsAsFactors = FALSE)
load(file.path(fpstate, 'states_albers.RData'))

# Convert all to Albers with "region" in the name.
aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
tnc <- spTransform(tnc, CRSobj = CRS(aea_crs))

tnc@data <- tnc@data %>%
  mutate(id = rownames(tnc@data), region = as.character(ECODE_NAME))

# Subset out the regions that are outside the US.
tnc_unique <- unique(allgeo$TNC)
tnc <- subset(tnc, region %in% tnc_unique)

# Clip TNC to US boundaries
goodusabounds <- gUnaryUnion(states_albers)
tncdat <- tnc@data
tnc <- gIntersection(tnc, goodusabounds, byid = TRUE, id = row.names(tnc@data))
tnc <- SpatialPolygonsDataFrame(tnc, tncdat)


# Loop through and make a plot, including both BBS and FIA locations, for the predictor variables @ 50 km
cols <- c(bg = 'gray80', text = 'black', state = 'gray20', fill = 'white')

# Draw the plots in a single figure with the 7 predictors
maps_geo <- allgeo_long %>%
  filter(!is.na(value)) %>%
  group_by(variable) %>%
  do(maps = point_map(., rbcol, tnc, states, bg_color=cols['bg'], text_color=cols['text'], state_color=cols['state'], fill_color = cols['fill'], show_legend = TRUE))

# Correct order of plots
geo_names_order <- c('temperature mean', 'precipitation mean', 'elevation diversity', 'soil diversity', 'geological age diversity', 'GPP diversity')

tw <- theme(plot.title = element_text(color = 'black'))
maps_geo$title <- geo_names_order

png(file.path(fpfig, 'geo_points.png'), height = 9, width = 9, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_geo$maps, maps_geo$title, function(p, name) ggplotGrob(p + ggtitle(name) + tw)), nrow = 3)
dev.off()

# Loop through and make a plot of the 9 response variables separately for BBS and for FIA

### bbs
maps_bbsbio <- bbsbio_long %>%
  filter(!is.na(value)) %>%
  group_by(variable) %>%
  do(maps = point_map(., rbcol, tnc, states, bg_color=cols['bg'], text_color=cols['text'], state_color=cols['state'], fill_color = cols['fill'], show_legend = TRUE))

tw <- theme(plot.title = element_text(color = 'black'))
maps_bbsbio$bio_title <- bio_titles[match(maps_bbsbio$variable, bio_names)] # Short title by bio variable
maps_bbsbio <- maps_bbsbio[match(bio_titles, maps_bbsbio$bio_title),] # Put bio variables in correct order
png(file.path(fpfig, 'bbs_response_points.png'), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_bbsbio$maps, maps_bbsbio$bio_title, function(p, name) ggplotGrob(p + ggtitle(name) + tw)), nrow = 3)
dev.off()

### fia
maps_fiabio <- fiabio_long %>%
  filter(!is.na(value)) %>%
  group_by(variable) %>%
  do(maps = point_map(., rbcol, tnc, states, bg_color=cols['bg'], text_color=cols['text'], state_color=cols['state'], fill_color = cols['fill'], show_legend = TRUE))

tw <- theme(plot.title = element_text(color = 'black'))
maps_fiabio$bio_title <- bio_titles[match(maps_fiabio$variable, bio_names)] # Short title by bio variable
maps_fiabio <- maps_fiabio[match(bio_titles, maps_fiabio$bio_title),] # Put bio variables in correct order
png(file.path(fpfig, 'fia_response_points.png'), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_fiabio$maps, maps_fiabio$bio_title, function(p, name) ggplotGrob(p + ggtitle(name) + tw)), nrow = 3)
dev.off()

