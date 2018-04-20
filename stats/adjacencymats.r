# Create adjacency matrix from HUC, BCR, and TNC ecoregion shape files.

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


# Bcr ---------------------------------------------------------------------

bcr <- subset(bcr, COUNTRY == 'USA')
library(rgeos)
bnames <- bcr@data$BCRNAME
bcr_combined <- gUnaryUnion(bcr, id = as.character(bnames))


bcr_neighb <- poly2nb(bcr_combined, row.names = names(bcr_combined), queen = TRUE)


# Tnc ---------------------------------------------------------------------

tnc <- subset(tnc, WWF_REALM2 == 'Nearctic' | ECO_NAME == 'Tropical Florida')
tnc_neighb <- poly2nb(tnc, row.names = tnc$ECODE_NAME, queen = TRUE)



# Make these into adj matrices --------------------------------------------

nb2binary <- function(neighb) {
  x <- matrix(0, nrow=length(neighb), ncol=length(neighb))
  for (i in 1:length(neighb)) {
    x[i, neighb[[i]]] <- 1
  }
  x
}

bcr_bin <- nb2binary(bcr_neighb)
dimnames(bcr_bin) <- list(names(bcr_combined), names(bcr_combined))

huc_bin <- nb2binary(huc_neighb)
dimnames(huc_bin) <- list(huc4$HUC4, huc4$HUC4)

tnc_bin <- nb2binary(tnc_neighb)
dimnames(tnc_bin) <- list(tnc$ECODE_NAME, tnc$ECODE_NAME)

save(bcr_bin, huc_bin, tnc_bin, file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/ecoregion_adjacency.RData')
