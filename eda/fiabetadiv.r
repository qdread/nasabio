# Test calculation of taxonomic beta diversity for FIA plots.

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'
fiapnw <- read.csv(file.path(fp, 'finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)

# Calculate basal area at plot level
library(dplyr)

fiasums_plot <- fiapnw %>%
  filter(STATUSCD == 1) %>%
  mutate(ba = pi * (DIA/200)^2) %>% # basal area in m2.
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT, MEASYEAR, SPCD) %>%
  summarize(basalarea = sum(ba),
            n = n())

fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

# Reproject FIA plots and calculate the distance in m from each plot to all other plots.
library(rgdal)

aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'

fiaspatial <- SpatialPointsDataFrame(coords = data.frame(x=fiacoords$lon, y=fiacoords$lat),
                                     data = as.data.frame(fiacoords[,1:4]),
                                     proj4string = CRS(wgs_crs)
                                    )

fiaalbers <- spTransform(fiaspatial, CRSobj = CRS(aea_crs))

# Fast function to get neighbor distances
# Refer to http://gis.stackexchange.com/questions/132384/distance-to-nearest-point-for-every-point-same-spatialpointsdataframe-in-r
getNeighbors <- function(dat, radius) {
  library(spdep)
  idlist <- dnearneigh(coordinates(dat), 0, radius)
  distlist <- nbdists(idlist, coordinates(dat))
  dflist <- list()
  for (i in 1:length(idlist)) {
    if (any(distlist[[i]] <= radius)) {
      dflist[[i]] <- data.frame(idx = idlist[[i]], dist = distlist[[i]][distlist[[i]] <= radius])
    }
    else {
      dflist[[i]] <- NA
    }
  }
  dflist
}

# Pairwise beta diversity between target point and everything within a chosen radius,

library(vegan)
library(vegetarian)

fianhb <- getNeighbors(fiaalbers, radius = 5000) # 5km radius
fia_shannonbetadiv5000 <- rep(NA,nrow(fiaalbers))
n_neighb5000 <- rep(0, nrow(fiaalbers))

pb <- txtProgressBar(0, nrow(fiaalbers), style=3)

for (i in 1:nrow(fiacoords)) {
  if (class(fianhb[[i]]) == 'data.frame') {
    # Subset out the data frame with the nearest neighbors
    plotcns <- fiaalbers[c(i, fianhb[[i]]$idx), ]$PLT_CN
    dat_i <- subset(fiasums_plot, PLT_CN %in% plotcns)
    # Convert into a site x species matrix
    sppids <- sort(unique(dat_i$SPCD))
    mat_i <- dat_i %>% group_by(PLT_CN) %>% do(x = sapply(sppids, function(z) sum(.$basalarea[.$SPCD == z])))
    mat_i <- do.call('rbind', mat_i$x)
    # Calculate beta-diversity for that matrix.
    if(!is.null(mat_i)) {
      fia_shannonbetadiv5000[i] <- d(abundances = mat_i, lev = 'beta', wts = FALSE, q = 1)
      n_neighb5000[i] <- nrow(mat_i) - 1
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# After successfully testing this for 5km radius, try it for a number of different radii.
radii <- c(1000,2000,3000,4000,5000,10000)
fia_shannonbetadiv <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_meanpairwisedissim <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_meanpairwisedissim_pa <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_phypairwise <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_phypairwise_pa <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_nneighb <- matrix(0, nrow = nrow(fiaalbers), ncol = length(radii))

pb2 <- txtProgressBar(0, length(radii) * nrow(fiaalbers), style = 3)
i <- 0

for (r in 1:length(radii)) {
  fianhb_r <- getNeighbors(fiaalbers, radius = radii[r])
  for (p in 1:nrow(fiaalbers)) {
    i <- i + 1
    if (class(fianhb_r[[p]]) == 'data.frame') {
      # Subset out the data frame with the nearest neighbors
      plotcns <- fiaalbers[c(p, fianhb_r[[p]]$idx), ]$PLT_CN
      dat_p <- subset(fiasums_plot, PLT_CN %in% plotcns)
      # Convert into a site x species matrix
      sppids <- sort(unique(dat_p$SPCD))
      mat_p <- dat_p %>% group_by(PLT_CN) %>% do(x = sapply(sppids, function(z) sum(.$basalarea[.$SPCD == z])))
      mat_p <- do.call('rbind', mat_p$x)
      # Calculate beta-diversity for that matrix.
      if(!is.null(mat_p)) {
        fia_shannonbetadiv[p, r] <- d(abundances = mat_p, lev = 'beta', wts = FALSE, q = 1)
        fia_meanpairwisedissim[p, r] <- mean(vegdist(x = mat_p, binary = FALSE, method = 'jaccard'))
        fia_meanpairwisedissim_pa[p, r] <- mean(vegdist(x = mat_p, binary = TRUE, method = 'jaccard'))
        fia_nneighb[p, r] <- nrow(mat_p) - 1
      }
    }
    setTxtProgressBar(pb2, i)
  }
}

close(pb2)

save(fia_shannonbetadiv, fia_nneighb, file = file.path(fp, 'fia_taxonomicbetadiv.RData'))
