# Create adjacency matrix from HUC, BCR, and TNC ecoregion shape files.
# Modified 25 Apr 2018: Get rid of Alaska and island ones that are not connected to others, and use built in function to produce binary matrix.
# Modified 08 Apr 2019: Create 10 "super-regions" for cross-validation from TNC ecoregions

fphuc <- '~/Dropbox/projects/aquaxterra/hucshapefiles'
fpregion <- '~/Dropbox/projects/nasabiodiv/regions'
fpfig <- '~/google_drive/NASABiodiversityWG/Figures/bbs_coefficient_maps'

library(rgdal)
library(dplyr)
library(maptools)
library(spdep)
huc4 <- readOGR(dsn = fphuc, layer = 'HU4_CONUS_Alb')
bcr <- readOGR(dsn = fpregion, layer = 'BCR_Terrestrial_master')
tnc <- readOGR(dsn = fpregion, layer = 'tnc_terr_ecoregions')


# Huc4 --------------------------------------------------------------------

huc4 <- spTransform(huc4, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
huc_neighb <- poly2nb(huc4, row.names = huc4$HUC4, queen = TRUE) # Long job.
huc_bin <- nb2mat(huc_neighb, style = 'B')

# Bcr ---------------------------------------------------------------------

bcr <- subset(bcr, COUNTRY == 'USA' & !PROVINCE_S %in% c('HAWAIIAN ISLANDS', 'ALASKA'))
library(rgeos)
bnames <- bcr@data$BCRNAME
bcr_combined <- gUnaryUnion(bcr, id = as.character(bnames))


bcr_neighb <- poly2nb(bcr_combined, row.names = names(bcr_combined), queen = TRUE)
bcr_bin <- nb2mat(bcr_neighb, style = 'B')

# Tnc ---------------------------------------------------------------------

tnc <- subset(tnc, (WWF_REALM2 == 'Nearctic' | ECO_NAME == 'Tropical Florida') & !grepl('Bermuda|Aleutian', ECODE_NAME))
tnc_neighb <- poly2nb(tnc, row.names = tnc$ECODE_NAME, queen = TRUE)
tnc_bin <- nb2mat(tnc_neighb, style = 'B')

save(bcr_bin, huc_bin, tnc_bin, file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/ecoregion_adjacency.RData')


# Get the ecoregion IDs where we actually have data -----------------------

load('~/Dropbox/projects/nasabiodiv/modelfits/bbs_spatial_mm_dat_50k.RData')
load('~/Dropbox/projects/nasabiodiv/modelfits/fia_spatial_mm_dat_50k.RData')
tnc_regions_use <- union(bbsgeo$TNC, fiageo$TNC)

# Create 10 super regions from TNC ----------------------------------------

# Use this code: https://github.com/BastienFR/clusteRneighbors
# See also: https://gis.stackexchange.com/questions/303524/clustering-spatially-connected-polygons-so-all-clusters-have-approximately-the-s

source('~/Dropbox/projects/nasabiodiv/code/bastien_rfunctions.r')
aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
fpregion <- '~/Dropbox/projects/nasabiodiv/regions'

library(sf)
library(dplyr)

tnc_st <- st_read(file.path(fpregion, 'tnc_terr_ecoregions.shp'))

# Get only the polygons that are in the US, calculate adjacency matrix and the centroid of each (in Albers)
tnc_st <- tnc_st %>% 
  st_transform(crs = aea_crs) %>%
  filter(ECODE_NAME %in% tnc_regions_use)
tnc_st_neigh <- st_touches(tnc_st, sparse = TRUE) %>% 
  setNames(1:nrow(.))
tnc_st_centroid <- st_centroid(tnc_st)

## make the first grouping with the required target size (it will be a minimum but for some exceptions see below))
set.seed(666)
gr1 <- cluster_neighbours(neighb_list = tnc_st_neigh, gr_size = 3, centro_pol = tnc_st_centroid)
sapply(gr1, length)

# Plot raw group
rawgrpdf <- data.frame(region = unlist(gr1), rawgrp = rep(1:length(gr1), unlist(sapply(gr1, length))))
tnc_st$rawgrp <- factor(rawgrpdf$rawgrp[order(rawgrpdf$region)])
plot(tnc_st['rawgrp'])

## Fuse the tiny groups to bigger groups, you can specify a new minimum smaller than the first one, 
##  higher this minimum is, bigger can be your final groups (potentially ~ gr_size + hard_min_size)
final_group <- fuse_tiny_group(initial_gr = gr1, init_sh = tnc_st, hard_min_size = 5)
sapply(final_group, length)
finalgrpdf <- data.frame(region = unlist(final_group), finalgrp = rep(1:length(final_group), unlist(sapply(final_group, length))))
tnc_st$finalgrp <- factor(finalgrpdf$finalgrp[order(finalgrpdf$region)])
plot(tnc_st['finalgrp'])

# See how many are in each.
final_group_names <- data.frame(TNC = tnc_st$ECODE_NAME, fold = factor(tnc_st$finalgrp)) # If using "fused" groups
#final_group_names <- data.frame(TNC = tnc_st$ECODE_NAME, fold = factor(tnc_st$rawgrp)) # If using the initial groups

bbsgeo %>% left_join(final_group_names) %>% group_by(fold) %>% summarize(n = n())
fiageo %>% left_join(final_group_names) %>% group_by(fold) %>% summarize(n = n())

write.csv(final_group_names, '~/Dropbox/projects/nasabiodiv/ecoregion_folds.csv', row.names = FALSE)
