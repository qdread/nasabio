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

# Make map with spplot instead of ggplot
model_map_spplot <- function(coefs, vars_to_plot, regions, state_borders, plot_title, ncut = 16) {
  require(latticeExtra)
  regions@data <- regions@data %>% 
    left_join(coefs)
  
  rbpal <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'RdYlBu')))(ncut+1)
  
  # Set custom cuts
  cuts_around_zero <- function(x,n) {
    xend <- max(abs(na.omit(x)))
    res <- seq(-xend, xend, length.out = n+1)
    res
  }
  
  state_layer <- list("sp.lines", state_borders, col = "gray50", first = FALSE)
  na_layer <- list("sp.polygons", regions, fill = "gray80", first = TRUE) # for NA values
  map_list <- list()
  for (i in vars_to_plot) {
    xcuts <- cuts_around_zero(regions@data[,i], ncut)
    
    map_list[[length(map_list) + 1]] <- 
      spplot(regions, zcol = i, col.regions = rbpal, at = xcuts, col = 'white',
             main = plot_title,
             sp.layout = list(state_layer, na_layer), 
             par.settings = list(panel.background=list(col="black")))
      
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
#load('C:/Users/Q/Dropbox/projects/nasabiodiv/mmfits.RData')

library(sp)
library(rgdal)
library(ggplot2)
library(dplyr)
library(purrr)
fpregion <- ('/mnt/research/nasabio/data/ecoregions')
#fpregion <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/regions'
huc4 <- readOGR(dsn = fpregion, layer = 'HU4_CONUS_Alb')
#huc4 <- readOGR(dsn = 'C:/Users/Q/Dropbox/projects/aquaxterra/hucshapefiles', layer = 'HU4_CONUS_Alb')
bcr <- readOGR(dsn = fpregion, layer = 'BCR_Terrestrial_master')
tnc <- readOGR(dsn = fpregion, layer = 'tnc_terr_ecoregions')
#load('~/states_albers.RData')
#load('~/R/states_albers.RData')
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

# Subset out the regions that are outside the US.
bcr <- subset(bcr, region %in% fit_bbs_bcr[[1]][[2]]$region & COUNTRY %in% 'USA' & !PROVINCE_S %in% 'ALASKA')
tnc <- subset(tnc, region %in% fit_bbs_tnc[[1]][[2]]$region)

# Create all maps ------------------------------------------

rbfill <- scale_fill_gradient2(low = "#4575B4", high = "#D73027", midpoint = 0)

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
prednames <- c('elevation_5k_100_sd', 'bio1_5k_100_mean', 'geological_age_5k_100_diversity', 'soil_type_5k_100_diversity', 'bio12_5k_100_mean', 'bio12_5k_100_sd', 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean')
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

### FIA
# Separate incidence and abundance
fpfig <- '/mnt/research/nasabio/figs/fia_coefficient_maps'

# Index the incidence based and abundance based variables
# Note that there is no beta pd and fd for FIA
fia_bio_names <- c("alpha_richness", "alpha_effspn", "alpha_phy_pa", "alpha_phy",
                   "alpha_func_pa", "alpha_func", "beta_td_sorensen_pa", "beta_td_sorensen",
                   "gamma_richness", "gamma_effspn", "gamma_phy_pa", "gamma_phy",
                   "gamma_func_pa", "gamma_func")
bio_titles <- c('alpha TD', 'alpha TD abundance', 'alpha PD', 'alpha PD abundance', 'alpha FD', 'alpha FD abundance', 'beta TD', 'beta TD abundance', 'gamma TD', 'gamma TD abundance', 'gamma PD', 'gamma PD abundance', 'gamma FD', 'gamma FD abundance')

bio_titles_incidence <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'gamma PD', 'alpha FD', 'gamma FD')
bio_titles_abundance <- paste(bio_titles_incidence, 'abundance')

maps_fia_huc <- map(fit_fia_huc, function(fit) model_map(fit$coef, prednames, rbfill, huc4, states))
maps_fia_bcr <- map(fit_fia_bcr, function(fit) model_map(fit$coef, prednames, rbfill, bcr, states))
maps_fia_tnc <- map(fit_fia_tnc, function(fit) model_map(fit$coef, prednames, rbfill, tnc, states))

maps_fia_huc_incidence <- maps_fia_huc[c(1,7,9,3,11,5,13)]
maps_fia_huc_abundance <- maps_fia_huc[c(2,8,10,4,12,6,14)]
maps_fia_bcr_incidence <- maps_fia_bcr[c(1,7,9,3,11,5,13)]
maps_fia_bcr_abundance <- maps_fia_bcr[c(2,8,10,4,12,6,14)]
maps_fia_tnc_incidence <- maps_fia_tnc[c(1,7,9,3,11,5,13)]
maps_fia_tnc_abundance <- maps_fia_tnc[c(2,8,10,4,12,6,14)]
map_mat <- rbind(c(1,2,3),
                 c(4,NA,5),
                 c(6,NA,7))

for (i in 1:8) {
  print(i)
  png(file.path(fpfig, paste0('HUC4_', geo_names[i], '_incidence.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_fia_huc_incidence, bio_titles_incidence, function(p, name) ggplotGrob(p[[i]] + ggtitle(name) + tw)), layout_matrix = map_mat)
  dev.off()
  png(file.path(fpfig, paste0('HUC4_', geo_names[i], '_abundance.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_fia_huc_abundance, bio_titles_abundance, function(p, name) ggplotGrob(p[[i]] + ggtitle(name) + tw)), layout_matrix = map_mat)
  dev.off()
  png(file.path(fpfig, paste0('BCR_', geo_names[i], '_incidence.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_fia_bcr_incidence, bio_titles_incidence, function(p, name) ggplotGrob(p[[i]] + ggtitle(name) + tw)), layout_matrix = map_mat)
  dev.off()
  png(file.path(fpfig, paste0('BCR_', geo_names[i], '_abundance.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_fia_bcr_abundance, bio_titles_abundance, function(p, name) ggplotGrob(p[[i]] + ggtitle(name) + tw)), layout_matrix = map_mat)
  dev.off()
  png(file.path(fpfig, paste0('TNC_', geo_names[i], '_incidence.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_fia_tnc_incidence, bio_titles_incidence, function(p, name) ggplotGrob(p[[i]] + ggtitle(name) + tw)), layout_matrix = map_mat)
  dev.off()
  png(file.path(fpfig, paste0('TNC_', geo_names[i], '_abundance.png')), height = 9, width = 12, res = 400, units = 'in')
  grid.arrange(grobs = map2(maps_fia_tnc_abundance, bio_titles_abundance, function(p, name) ggplotGrob(p[[i]] + ggtitle(name) + tw)), layout_matrix = map_mat)
  dev.off()
}

# Create all maps (sp version) --------------------------------------------

### BBS
prednames <- c('elevation_5k_100_sd', 'bio1_5k_100_mean', 'geological_age_5k_100_diversity', 'soil_type_5k_100_diversity', 'bio12_5k_100_mean', 'bio12_5k_100_sd', 'dhi_gpp_5k_100_sd', 'human_footprint_5k_100_mean')
fpfig <- '/mnt/research/nasabio/figs/bbs_coefficient_maps'
#fpfig <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/figs/bbs_coefficient_maps'
bio_titles <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'beta PD', 'gamma PD', 'alpha FD', 'beta FD', 'gamma FD')
geo_names <- c('elevation_sd','temperature_mean','geol_age_diversity','soil_diversity','precip_mean','precip_sd','gpp_sd','footprint_mean')

# Create nested list of maps for each predictor by response combo (9 predictors x 8 responses)
maps_bbs_huc <- map2(fit_bbs_huc, bio_titles, function(fit, name) model_map_spplot(fit$coef, prednames, huc4, states_albers, name))
maps_bbs_bcr <- map2(fit_bbs_bcr, bio_titles, function(fit, name) model_map_spplot(fit$coef, prednames, bcr, states_albers, name))
maps_bbs_tnc <- map2(fit_bbs_tnc, bio_titles, function(fit, name) model_map_spplot(fit$coef, prednames, tnc, states_albers, name))

# Create grids of maps by predictor, and grids of maps by response
library(gridExtra)

for (i in 1:8) {
  print(i)
  png(file.path(fpfig, paste0('HUC4_', geo_names[i], '.png')), height = 9, width = 12, res = 400, units = 'in')
    grid.arrange(grobs = map(maps_bbs_huc, i), nrow = 3)
  dev.off()
  png(file.path(fpfig, paste0('BCR_', geo_names[i], '.png')), height = 9, width = 12, res = 400, units = 'in')
    grid.arrange(grobs = map(maps_bbs_bcr, i), nrow = 3)
  dev.off()
  png(file.path(fpfig, paste0('TNC_', geo_names[i], '.png')), height = 9, width = 12, res = 400, units = 'in')
    grid.arrange(grobs = map(maps_bbs_tnc, i), nrow = 3)
  dev.off()
}

### FIA
# Separate incidence and abundance
fpfig <- '/mnt/research/nasabio/figs/fia_coefficient_maps'

# Index the incidence based and abundance based variables
# Note that there is no beta pd and fd for FIA
fia_bio_names <- c("alpha_richness", "alpha_effspn", "alpha_phy_pa", "alpha_phy",
                   "alpha_func_pa", "alpha_func", "beta_td_sorensen_pa", "beta_td_sorensen",
                   "gamma_richness", "gamma_effspn", "gamma_phy_pa", "gamma_phy",
                   "gamma_func_pa", "gamma_func")
bio_titles <- c('alpha TD', 'alpha TD abundance', 'alpha PD', 'alpha PD abundance', 'alpha FD', 'alpha FD abundance', 'beta TD', 'beta TD abundance', 'gamma TD', 'gamma TD abundance', 'gamma PD', 'gamma PD abundance', 'gamma FD', 'gamma FD abundance')

maps_fia_huc <- map2(fit_fia_huc, bio_titles, function(fit, name) model_map_spplot(fit$coef, prednames, huc4, states_albers, name))
maps_fia_bcr <- map2(fit_fia_bcr, bio_titles, function(fit, name) model_map_spplot(fit$coef, prednames, bcr, states_albers, name))
maps_fia_tnc <- map2(fit_fia_tnc, bio_titles, function(fit, name) model_map_spplot(fit$coef, prednames, tnc, states_albers, name))

maps_fia_huc_incidence <- maps_fia_huc[c(1,7,9,3,11,5,13)]
maps_fia_huc_abundance <- maps_fia_huc[c(2,8,10,4,12,6,14)]
maps_fia_bcr_incidence <- maps_fia_bcr[c(1,7,9,3,11,5,13)]
maps_fia_bcr_abundance <- maps_fia_bcr[c(2,8,10,4,12,6,14)]
maps_fia_tnc_incidence <- maps_fia_tnc[c(1,7,9,3,11,5,13)]
maps_fia_tnc_abundance <- maps_fia_tnc[c(2,8,10,4,12,6,14)]
map_mat <- rbind(c(1,2,3),
                 c(4,NA,5),
                 c(6,NA,7))

for (i in 1:8) {
  print(i)
  png(file.path(fpfig, paste0('HUC4_', geo_names[i], '_incidence.png')), height = 9, width = 12, res = 400, units = 'in')
	grid.arrange(grobs = map(maps_fia_huc_incidence, i), layout_matrix = map_mat)
  dev.off()
  png(file.path(fpfig, paste0('BCR_', geo_names[i], '_incidence.png')), height = 9, width = 12, res = 400, units = 'in')
	grid.arrange(grobs = map(maps_fia_huc_abundance, i), layout_matrix = map_mat)
  dev.off()
  png(file.path(fpfig, paste0('TNC_', geo_names[i], '_incidence.png')), height = 9, width = 12, res = 400, units = 'in')
	grid.arrange(grobs = map(maps_fia_bcr_incidence, i), layout_matrix = map_mat)
  dev.off()
  png(file.path(fpfig, paste0('HUC4_', geo_names[i], '_abundance.png')), height = 9, width = 12, res = 400, units = 'in')
	grid.arrange(grobs = map(maps_fia_bcr_abundance, i), layout_matrix = map_mat)
  dev.off()
  png(file.path(fpfig, paste0('BCR_', geo_names[i], '_abundance.png')), height = 9, width = 12, res = 400, units = 'in')
	grid.arrange(grobs = map(maps_fia_bcr_incidence, i), layout_matrix = map_mat)
  dev.off()
  png(file.path(fpfig, paste0('TNC_', geo_names[i], '_abundance.png')), height = 9, width = 12, res = 400, units = 'in')
	grid.arrange(grobs = map(maps_fia_bcr_abundance, i), layout_matrix = map_mat)
  dev.off()
}