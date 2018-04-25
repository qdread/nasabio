# Create adjacency matrix from HUC, BCR, and TNC ecoregion shape files.
# Modified 25 Apr 2018: Get rid of Alaska and island ones that are not connected to others, and use built in function to produce binary matrix.


fphuc <- 'C:/Users/Q/Dropbox/projects/aquaxterra/hucshapefiles'
fpregion <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/regions'
fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_coefficient_maps'

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
