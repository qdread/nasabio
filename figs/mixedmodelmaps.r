# Create shaded polygon maps for mixed model fits
# QDR NASABIOXGEO 08 Mar 2018


# Define functions --------------------------------------------------------

# Make map (run fortify each time this function is called, slower but simpler)
# Return the plots as a list so that we can tile them in different ways
model_map <- function(coefs, vars_to_plot, fill_scale, regions, state_borders) {
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
          legend.title = element_blank())
  
  map_list <- list()
  for (i in vars_to_plot) {
    map_list[[length(map_list) + 1]] <- ggplot(region_fort) +
      geom_polygon(aes_string(x='long', y='lat', group='group', fill=i)) +
      geom_path(aes_string(x='long', y='lat', group='group'), color = 'white', size = 0.25) +
      geom_path(data = state_borders, aes(x = long, y = lat, group = group), color = 'gray20') +
      fill_scale +
      coord_equal() +
      blktheme
    print(i)
  }
  return(map_list)
}

# Make predicted vs observed plot
obs_pred_plot <- function() {
  
}

# Load region data --------------------------------------------------------

# Load fits
load('/mnt/research/nasabio/temp/mmfits.RData')

library(rgdal)
library(ggplot2)
fpregion <- file.path(fp, 'ecoregions')
huc4 <- readOGR(dsn = fpregion, layer = 'HU4_CONUS_Alb')
bcr <- readOGR(dsn = fpregion, layer = 'BCR_Terrestrial_master')
tnc <- readOGR(dsn = fpregion, layer = 'tnc_terr_ecoregions')
states <- read.csv('~/states_albers.csv', stringsAsFactors = FALSE)

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


# Create all maps ------------------------------------------

rbfill <- scale_fill_gradient2(low = "#4575B4", high = "#D73027", midpoint = 0)

# Subset out the regions that are outside the US.
bcr <- subset(bcr, region %in% fit_bbs_bcr[[1]][[2]]$region & COUNTRY %in% 'USA' & !PROVINCE_S %in% 'ALASKA')
tnc <- subset(tnc, region %in% fit_bbs_tnc[[1]][[2]]$region)

# Create nested list of maps for each predictor by response combo (9 predictors x 8 responses)
maps_bbs_huc <- map(fit_bbs_huc, function(fit) model_map(fit$coef, prednames, rbfill, huc4, states))
maps_bbs_bcr <- map(fit_bbs_bcr, function(fit) model_map(fit$coef, prednames, rbfill, bcr, states))
maps_bbs_tnc <- map(fit_bbs_tnc, function(fit) model_map(fit$coef, prednames, rbfill, tnc, states))

names(maps_bbs_huc) <- names(bbsbio)[-1]
maps_bbs_huc <- lapply(maps_bbs_huc, setNames, nm = prednames)
names(maps_bbs_bcr) <- names(bbsbio)[-1]
maps_bbs_bcr <- lapply(maps_bbs_bcr, setNames, nm = prednames)
names(maps_bbs_tnc) <- names(bbsbio)[-1]
maps_bbs_tnc <- lapply(maps_bbs_tnc, setNames, nm = prednames)

# Create grids of maps by predictor, and grids of maps by response
library(gridExtra)
fpfig <- '/mnt/research/nasabio/figs/bbs_coefficient_maps'
bio_titles <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'beta PD', 'gamma PD', 'alpha FD', 'beta FD', 'gamma FD')
geo_names <- c('elevation_sd','temperature_mean','geol_age_diversity','soil_diversity','precip_mean','precip_sd','gpp_sd','footprint_mean')
tw <- theme(plot.title = element_text(color = 'white'))
for (i in 1:8) {
  print(i)
  png(file.path(fpfig, paste0('HUC4_', geo_names[i], '.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_bbs_huc, bio_titles, function(p, name) ggplotGrob(p[[i]] + ggtitle(name) + tw)), nrow = 3)
  dev.off()
  png(file.path(fpfig, paste0('BCR_', geo_names[i], '.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_bbs_bcr, bio_titles, function(p, name) ggplotGrob(p[[i]] + ggtitle(name) + tw)), nrow = 3)
  dev.off()
  png(file.path(fpfig, paste0('TNC_', geo_names[i], '.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_bbs_tnc, bio_titles, function(p, name) ggplotGrob(p[[i]] + ggtitle(name) + tw)), nrow = 3)
  dev.off()
}

