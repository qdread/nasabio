# Create shaded polygon maps for mixed model fits
# With black background, for ESA 2019 PPT
# QDR NASABIOXGEO 18 June 2019 (convert to run locally)

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
      geom_polygon(data = region_fort %>% filter(grepl('Black Hills', region)), aes(x=long, y=lat, group=group, fill=Estimate)) +
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
fpcoef <- '~/Dropbox/projects/nasabiodiv/modelfits' # Local
fpregion <- '~/Dropbox/projects/nasabiodiv/regions'
fpstate <- '~/Documents/R'


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
tnc_unique <- na.omit(unique(coef_all$region[coef_all$ecoregion == 'TNC']))
tnc <- subset(tnc, region %in% tnc_unique)

# Clip TNC to US boundaries
goodusabounds <- gUnaryUnion(states_albers)
tncdat <- tnc@data
tnc <- gIntersection(tnc, goodusabounds, byid = TRUE, id = tnc$id)
row.names(tncdat) <- tncdat$id
tnc <- SpatialPolygonsDataFrame(tnc, tncdat)

# Create all maps ------------------------------------------

# function to produce a fill scale with the maximum and minimum evenly divided above and below zero.
rbfill <- function(x) {
	xlimit <- max(abs(x))
	xlimit_roundup <- signif(xlimit + 0.1, 2)
	xbreaks <- c(-xlimit_roundup, 0, xlimit_roundup)
	scale_fill_gradient2(low = "#4575B4", high = "#D73027", midpoint = 0, limits = c(-xlimit_roundup, xlimit_roundup), labels = xbreaks, breaks = xbreaks)
}
signif_fill <- scale_fill_manual(values = c(negative = "#4575B4", positive = "#D73027", zero = 'gray50'))

library(gridExtra)
fpfig <- '~/google_drive/NASABiodiversityWG/Conferences/ESA2019/talkimgs/maps'

bio_titles <- c('alpha taxonomic', 'alpha phylogenetic', 'alpha functional', 'beta taxonomic', 'beta phylogenetic', 'beta functional', 'gamma taxonomic', 'gamma phylogenetic', 'gamma functional')
bio_names <- c("alpha_richness", "alpha_phy_pa", "alpha_func_pa",
               "beta_td_sorensen_pa", "beta_phy_pa", "beta_func_pa",
               "gamma_richness", "gamma_phy_pa", "gamma_func_pa")
geo_names <- c('elevation_diversity','temperature_mean','geol_age_diversity','soil_diversity','precip_mean','gpp_sd', 'intercept')
prednames <- c('elevation_5k_tri_50_mean', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'dhi_gpp_5k_tri_50_mean', 'Intercept')

cols <- c(bg = 'black', text = 'white', state = 'gray20')

# Just make some taxonomic diversity maps

  maps_bbs_coef <- coef_wide %>%
    filter(taxon == 'bbs', response %in% c('alpha_richness', 'beta_td_sorensen_pa', 'gamma_richness'), parameter %in% c('elevation_5k_tri_50_mean', 'dhi_gpp_5k_tri_50_mean', 'bio12_5k_50_mean')) %>%
    dplyr::select(response, parameter, region, Estimate) %>%
    rename(ECODE_NAME = region) %>%
    group_by(response, parameter) %>%
    do(maps = model_map(., rbfill(.$Estimate), tnc, states, bg_color=cols['bg'], text_color=cols['text'], state_color=cols['state']))
  maps_fia_coef <- coef_wide %>%
    filter(taxon == 'fia',  response %in% c('alpha_richness', 'beta_td_sorensen_pa', 'gamma_richness'), parameter %in% c('elevation_5k_tri_50_mean', 'dhi_gpp_5k_tri_50_mean', 'bio12_5k_50_mean')) %>%
    dplyr::select(response, parameter, region, Estimate) %>%
    rename(ECODE_NAME = region) %>%
    group_by(response, parameter) %>%
    do(maps = model_map(., rbfill(.$Estimate), tnc, states, bg_color=cols['bg'], text_color=cols['text'], state_color=cols['state']))	

# Write them all to separate files with no text in the image.
pwalk(maps_bbs_coef, function(response, parameter, maps) ggsave(file.path(fpfig, paste0('bbs_', response, '_', parameter, '.png')), maps, height = 3, width = 4, dpi = 400))
pwalk(maps_fia_coef, function(response, parameter, maps) ggsave(file.path(fpfig, paste0('fia_', response, '_', parameter, '.png')), maps, height = 3, width = 4, dpi = 400))


# Maps with point biodiv and geodiv values --------------------------------

# Modified from maps put in the 2018 iale presentation
# Change: get rid of text in legend, make legend numbers bigger
library(sp)
library(rgdal)
library(ggplot2)
library(dplyr)
library(purrr)
fp <- '~/Dropbox/projects/nasabiodiv/'
fpfig <- '~/google_drive/NASABiodiversityWG/Conferences/ESA2019/talkimgs/maps'
states <- read.csv('~/Documents/R/states_albers.csv', stringsAsFactors = FALSE)
bbscoords <- read.csv(file.path(fp,'bbs_correct_route_centroids.csv'))

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

# BBS
bbsbio <- read.csv(file.path(fp, 'bbs_allbio_wide.csv'), stringsAsFactors = FALSE)
bbsgeo <- read.csv(file.path(fp, 'bbs_allgeo_wide.csv'), stringsAsFactors = FALSE)

bbsgeomeans <- bbsgeo %>%
  dplyr::select(rteNo, elevation_5k_50_mean, geological_age_5k_50_mode, dhi_gpp_5k_50_mean)

bbsbio <- bbsbio %>%
  dplyr::select(rteNo, alpha_richness_50, beta_td_sorensen_pa_50, gamma_richness_50)
bbsgeo <- bbsgeo %>%
  dplyr::select(rteNo, elevation_5k_tri_50_mean, geological_age_5k_50_diversity, dhi_gpp_5k_tri_50_mean)

bbscoords <- bbscoords %>% left_join(bbsbio) %>% left_join(bbsgeo)
cscale <- function(b, n) scale_color_gradientn(name = n, colours = colorRampPalette(RColorBrewer::brewer.pal(9,'YlOrRd'), bias = b)(50))
biases <- c(1, 2, 0.8, 2, 1, 1)
legnames <- c('Alpha diversity', 'Beta diversity', 'Gamma diversity', 'Elevation variability', 'Geological age diversity', 'Productivity variability')

for (i in 6:8) {
  dati <- bbscoords %>% 
    filter(lat < 50, !is.na(bbscoords[,i]))
  p <- ggplot() +
    #geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.75) +
    geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'white', size = 0.75) +
    geom_point(data = dati, aes_string(x='lon.1', y='lat.1', color = names(bbscoords)[i]), size = 1) +
    coord_equal() +
    blktheme + cscale(biases[i-5], legnames[i-5]) + 
    theme(legend.text = element_text(color='white'),
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black', fill = 'black'),
          legend.position = c(0.15,0.1), 
          legend.key.width = unit(0.2, 'inches'),
          legend.direction = 'horizontal')
  ggsave(file.path(fpfig, paste0('point_bbs', names(bbscoords)[i], 'map.png')), p, height = 3, width = 4, dpi = 400)
}

# FIA
fiacoords <- read.csv(file.path(fp, 'fia_unfuzzed/fia_fuzzed_coords.csv'))

load(file.path(fp, 'modelfits/fia_spatial_mm_dat_50k.RData'))

fiageomeans <- fiageo %>%
    dplyr::select(PLT_CN, elevation_5k_50_mean,  dhi_gpp_5k_50_mean)

fiabio <- fiabio %>%
  dplyr::select(PLT_CN, alpha_richness, beta_td_sorensen_pa, gamma_richness)
fiageo <- fiageo %>%
  dplyr::select(PLT_CN, elevation_5k_tri_50_mean, geological_age_5k_50_diversity, dhi_gpp_5k_tri_50_mean)

# Combine bbs and fia geo
fiacoordsbio <- fiacoords %>% left_join(fiabio) %>% filter(!is.na(alpha_richness))
allgeocoords <- fiacoords %>% left_join(fiageo) %>% bind_rows(bbscoords) %>% filter(!is.na(elevation_5k_tri_50_mean))

# Combine bbs and gia geo means
bbsgeomeans <- bbscoords %>% left_join(bbsgeomeans)
allgeomeans <- fiacoords %>% left_join(fiageomeans) %>% bind_rows(bbsgeomeans) %>% filter(!is.na(elevation_5k_50_mean))

biases2 <- c(1, 1, 0.8, 2, 1, 1)

for (i in 6:8) {
  dati <- fiacoordsbio %>% 
    filter(lat_fuzzed < 50, !is.na(fiacoordsbio[,i]))
  p <- ggplot() +
    #geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.75) +
    geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'white', size = 0.75) +
    geom_point(data = dati, aes_string(x='lon_aea_fuzzed', y='lat_aea_fuzzed', color = names(fiacoordsbio)[i]), size = 1) +
    coord_equal() +
    blktheme + cscale(biases2[i-5], legnames[i-5]) + 
    theme(legend.text = element_text(color='white'),
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black', fill = 'black'),
          legend.position = c(0.15,0.1), 
          legend.key.width = unit(0.2, 'inches'),
          legend.direction = 'horizontal')
  ggsave(file.path(fpfig, paste0('point_fia', names(fiacoordsbio)[i], 'map.png')), p, height = 3, width = 4, dpi = 400)
}

for (i in 6:8) {
  dati <- allgeocoords %>% 
    filter(lat < 50 | lat_fuzzed < 50, !is.na(allgeocoords[,i])) %>%
    mutate(lonaea = pmin(lon_aea_fuzzed, lon.1, na.rm=T),
           lataea = pmin(lat_aea_fuzzed, lat.1, na.rm=T))
  p <- ggplot() +
    #geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.75) +
    geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'white', size = 0.75) +
    geom_point(data = dati, aes_string(x='lonaea', y='lataea', color = names(allgeocoords)[i]), size = 1) +
    coord_equal() +
    blktheme + cscale(biases[i-2], legnames[i-2]) + 
    theme(legend.text = element_text(color='white'),
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black', fill = 'black'),
          legend.position = c(0.15,0.1), 
          legend.key.width = unit(0.2, 'inches'),
          legend.direction = 'horizontal')
  ggsave(file.path(fpfig, paste0('point_geo', names(allgeocoords)[i], 'map.png')), p, height = 3, width = 4, dpi = 400)
}

brks <- list(c(1000,2000,3000), c(10000, 20000))
cscale2 <- function(b, n, brks) scale_color_gradientn(name = n, breaks = brks, colours = colorRampPalette(RColorBrewer::brewer.pal(9,'YlOrRd'), bias = b)(50))

 # Maps with means
for (i in c(6,7)) {
  dati <- allgeomeans %>% 
    filter(lat < 50 | lat_fuzzed < 50, !is.na(allgeomeans[,i])) %>%
    mutate(lonaea = pmin(lon_aea_fuzzed, lon.1, na.rm=T),
           lataea = pmin(lat_aea_fuzzed, lat.1, na.rm=T))
  p <- ggplot() +
    #geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.75) +
    geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'white', size = 0.75) +
    geom_point(data = dati, aes_string(x='lonaea', y='lataea', color = names(allgeomeans)[i]), size = 1) +
    coord_equal() +
    blktheme + cscale2(1, legnames[i-2], brks[[i-5]]) + 
    theme(legend.text = element_text(color='white'),
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black', fill = 'black'),
          legend.position = c(0.15,0.1), 
          legend.key.width = unit(0.2, 'inches'),
          legend.direction = 'horizontal')
  ggsave(file.path(fpfig, paste0('pointmean_geo', names(allgeomeans)[i], 'map.png')), p, height = 3, width = 4, dpi = 400)
}

# Mode of geological age
dati <- allgeomeans %>% 
  filter(lat < 50 | lat_fuzzed < 50, !is.na(geological_age_5k_50_mode)) %>%
  mutate(lonaea = pmin(lon_aea_fuzzed, lon.1, na.rm=T),
         lataea = pmin(lat_aea_fuzzed, lat.1, na.rm=T))
p <- ggplot() +
  #geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.75) +
  geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'white', size = 0.75) +
  geom_point(data = dati, aes(x=lonaea, y=lataea, color = factor(geological_age_5k_50_mode)), size = 1) +
  coord_equal() +
  blktheme + 
  theme(legend.text = element_text(color='white'),
        legend.title = element_blank(),
        legend.background = element_rect(color = 'black', fill = 'black'),
        legend.position = 'none', 
        legend.key.width = unit(0.2, 'inches'),
        legend.direction = 'horizontal')
ggsave(file.path(fpfig, 'modegeoage.png'), p, height = 3, width = 4, dpi = 400)
