# Compile FIA environmental and tree data from California for Andrew
# QDR 07 May 2018

# Get subset of plot IDs in California

fiadat <- read.csv('/mnt/research/nasabio/data/fia/treedata10nov/finley_trees_continental_US_most_recent_evaluations_nov8_2017.csv')
fia_cali <- subset(fiadat, STATECD==6)

# Create a table by species code and basal area

library(dplyr)

fia_cali_basalarea <- fia_cali %>%
  filter(STATUSCD == 1) %>%
  mutate(ba = pi * (DIA/200)^2) %>% # basal area in m2.
  group_by(COUNTYCD, PLT_CN, PLOT, MEASYEAR, SPCD) %>%
  summarize(basalarea = sum(ba),
            n = n())
			
# Make a wide format table by basal area and stem number
library(reshape2)

fia_basal_wide <- dcast(fia_cali_basalarea, COUNTYCD + PLT_CN + PLOT + MEASYEAR ~ SPCD, value.var = 'basalarea')
fia_stemnumber_wide <- dcast(fia_cali_basalarea, COUNTYCD + PLT_CN + PLOT + MEASYEAR ~ SPCD, value.var = 'n')

# Get the true coordinates of the FIA plots and extract all bands from sierra_srad.tif at those locations

fiacoords <- read.csv('~/data/allfia.csv')
cali_plots <- unique(fia_cali$PLT_CN)
cali_coords <- subset(fiacoords, CN %in% cali_plots)

library(purrr)
library(raster)

# Get projection of the raster and reproject the plots into that projection
raster_file <- '/mnt/research/nasabio/data/climate/sierra_srad.tif'
ras <- brick(raster_file)
raster_crs <- proj4string(ras)

cali_sp <- SpatialPointsDataFrame(coords = cali_coords[,c('ACTUAL_LON', 'ACTUAL_LAT')], data = cali_coords[,'CN',drop=FALSE], proj4string = CRS('+proj=longlat +ellps=WGS84 +no_defs'))
cali_albers <- spTransform(cali_sp, CRSobj = CRS(raster_crs))

# Try this with the extract() function
insol_values_extract <- extract(ras, cali_albers) # Works better
insol_values_extract <- cbind(PLT_CN = cali_coords$CN, insol_values_extract)
insol_values_extract <- as.data.frame(insol_values_extract)

# Get whatever other point values I already extracted
fia_bypoint <- read.csv('/mnt/research/nasabio/data/fia/fia_geo_by_point.csv')
cali_bypoint <- subset(fia_bypoint, PLT_CN %in% cali_plots)

# Write csvs
fp <- '/mnt/research/nasabio/data/fia/sierra'
write.csv(fia_basal_wide, file = file.path(fp, 'cali_basalarea.csv'), row.names = FALSE)
write.csv(fia_stemnumber_wide, file = file.path(fp, 'cali_stemnumber.csv'), row.names = FALSE)
write.csv(cali_bypoint, file = file.path(fp, 'cali_environmentaldata.csv'), row.names = FALSE)
write.csv(insol_values_extract, file = file.path(fp, 'sierra_rad.csv'), row.names = FALSE)