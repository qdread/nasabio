# Make maps of region boundaries just to see what they look like

fphuc <- 'C:/Users/Q/Dropbox/projects/aquaxterra/hucshapefiles'
fpregion <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/regions'
fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_coefficient_maps'

library(rgdal)
library(dplyr)
library(ggplot2)
library(MapColoring) # To create four color fill to make the map look nice.
huc4 <- readOGR(dsn = fphuc, layer = 'HU4_CONUS_Alb')
bcr <- readOGR(dsn = fpregion, layer = 'BCR_Terrestrial_master')
tnc <- readOGR(dsn = fpregion, layer = 'tnc_terr_ecoregions')

# four color pastel-ish from brewer
cs <- scale_fill_manual(values = c('#7fc97f','#beaed4','#fdc086','#ffff99','#bf5b17','#f0027f'), guide = FALSE)

# Huc4 --------------------------------------------------------------------


huc4 <- spTransform(huc4, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
huc4colors <- getColoring(huc4)

huc4@data <- huc4@data %>% 
  mutate(id = rownames(huc4@data), HUC4 = as.character(HUC4), mapcolor = factor(huc4colors)) 

huc4_fort <- fortify(huc4, region = 'id') %>% left_join(huc4@data, by = 'id')

huc_map <- ggplot(huc4_fort) +
  geom_polygon(aes_string(x='long', y='lat', group='group', fill='mapcolor')) +
  geom_path(aes_string(x='long', y='lat', group='group'), color = 'white', size = 0.25) +
  borders('state', colour = 'gray20') +
  coord_map(xlim = c(-125,-67), ylim = c(25,50)) +
  theme_bw() + cs +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'black'), panel.border = element_blank(), plot.background = element_rect(fill = 'black'))
ggsave(file.path(fpfig, 'borders_huc4.png'), huc_map, height = 6, width = 9, dpi = 400)  


# BCR ---------------------------------------------------------------------

# Dissolve BCR regions that are split by state and country into single regions
# First get rid of regions outside the USA then run gUnaryUnion on the result.

bcr <- subset(bcr, COUNTRY == 'USA')
library(rgeos)
bnames <- bcr@data$BCRNAME
bcr_combined <- gUnaryUnion(bcr, id = as.character(bnames))

bcrcolors <- getColoring(bcr_combined)
bcrcolors <- data.frame(mapcolor = factor(bcrcolors), BCR = names(bcr_combined))

bcr@data <- bcr@data %>% 
  mutate(id = rownames(bcr@data), BCR = as.character(BCRNAME)) %>% left_join(bcrcolors)

bcr_fort <- fortify(bcr, region = 'id') %>% left_join(bcr@data, by = 'id')

bcr_map <- ggplot(bcr_fort) +
  geom_polygon(aes_string(x='long', y='lat', group='group', fill='mapcolor')) +
  geom_path(aes_string(x='long', y='lat', group='group'), color = 'white', size = 0.25) +
  borders('state', colour = 'gray20') +
  coord_map(xlim = c(-125,-67), ylim = c(25,50)) +
  theme_bw() + cs +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'black'), panel.border = element_blank(), plot.background = element_rect(fill = 'black'))
ggsave(file.path(fpfig, 'borders_bcr.png'), bcr_map, height = 6, width = 9, dpi = 400)


# TNC ---------------------------------------------------------------------

tnc <- subset(tnc, WWF_REALM2 == 'Nearctic' | ECO_NAME == 'Tropical Florida')

tnccolors <- getColoring(tnc)
tnc@data <- tnc@data %>% 
  mutate(id = rownames(tnc@data), TNC = as.character(ECODE_NAME), mapcolor = factor(tnccolors))

tnc_fort <- fortify(tnc, region = 'id') %>% left_join(tnc@data, by = 'id')

tnc_map <- ggplot(tnc_fort) +
  geom_polygon(aes_string(x='long', y='lat', group='group', fill='mapcolor')) +
  geom_path(aes_string(x='long', y='lat', group='group'), color = 'white', size = 0.25) +
  borders('state', colour = 'gray20') +
  coord_map(xlim = c(-125,-67), ylim = c(25,50)) +
  theme_bw() + cs +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'black'), panel.border = element_blank(), plot.background = element_rect(fill = 'black'))
ggsave(file.path(fpfig, 'borders_tnc.png'), tnc_map, height = 6, width = 9, dpi = 400)