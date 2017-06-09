# FIA taxonomic and phylogenetic beta diversity calculation for the cluster.
# Forked 17 April. Just preparatory work, actual diversity calculated in a different script.
# Modified 17 May: add rows to matrix that have zero species.

fp <- '/mnt/research/nasabio'
fiapnw <- read.csv(file.path(fp, 'data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)

#phylogenetic distance matrix
library(ape)
load(file.path(fp, 'data/fia/pnwphylo_potter.r'))
fiadist <- cophenetic(pnwphylo)

# Get FIA species code to scientific name lookup table
fiataxa <- read.csv(file.path(fp, 'data/fia/fia_taxon_lookuptable.csv'), stringsAsFactors = FALSE)
pnw_codes <- unique(fiapnw$SPCD)

# Convert fiapnw into a site x species matrix at plot level

library(dplyr)

area_by_sp <- function(dat, sppids) {
  areas <- numeric(length(sppids))
  for (i in 1:length(sppids)) {
    areas[i] <- sum(dat$basalarea[dat$SPCD == sppids[i]])
  }
  areas
}

# Plot level
fiasums_plot <- fiapnw %>%
  filter(STATUSCD == 1) %>%
  mutate(ba = pi * (DIA/200)^2) %>% # basal area in m2.
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT, MEASYEAR, SPCD) %>%
  summarize(basalarea = sum(ba),
            n = n())

sppids <- sort(unique(fiasums_plot$SPCD))

idx <- match(sppids, fiataxa$FIA.Code)
fiataxa$sciname <- paste(gsub(' ', '', fiataxa$Genus), gsub(' ', '', fiataxa$Species), sep = '_')
fiataxa$sciname[fiataxa$sciname == 'Chrysolepis_chrysophyllavar.chrysophylla'] <- 'Chrysolepis_chrysophylla'
fiataxa$sciname[fiataxa$sciname == 'Populus_balsamiferassp.Trichocarpa'] <- 'Populus_balsamifera_trichocarpa'

fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

#fiaplotlist <- fiasums_plot %>% do(x = area_by_sp(., sppids))
#fiaplotmat <- do.call('rbind', fiaplotlist$x)

fiaplotlist <- list()
pb <- txtProgressBar(0, nrow(fiacoords), style=3)
for (i in 1:nrow(fiacoords)) {
	fiaplotlist[[i]] <- area_by_sp(subset(fiasums_plot, PLT_CN == fiacoords$PLT_CN[i]), sppids)
	setTxtProgressBar(pb,i)
}
close(pb)
fiaplotmat <- do.call('rbind', fiaplotlist)

# Get the species names from the lookup table that go with the numerical codes.
#sppnames <- pnw_species$sciname[match(sppids, pnw_codes)]
dimnames(fiaplotmat)[[2]] <- fiataxa$sciname[idx]

# Add Abies shastensis to Abies magnifica
fiaplotmat[,'Abies_magnifica'] <- fiaplotmat[,'Abies_magnifica'] + fiaplotmat[,'Abies_shastensis']
fiaplotmat <- fiaplotmat[, !grepl('Abies_shastensis', dimnames(fiaplotmat)[[2]])]

# Get rid of the unknown species.
fiaplotmat <- fiaplotmat[, dimnames(fiaplotmat)[[2]] %in% pnwphylo$tip.label]


# Reproject FIA plots and calculate the distance in m from each plot to all other plots.
library(rgdal)

aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'

fiaspatial <- SpatialPointsDataFrame(coords = data.frame(x=fiacoords$lon, y=fiacoords$lat),
                                     data = as.data.frame(fiacoords[,1:4]),
                                     proj4string = CRS(wgs_crs)
)

fiaalbers <- spTransform(fiaspatial, CRSobj = CRS(aea_crs))

# Generate distance matrix for functional beta-diversity
source('~/code/fia/trydistmat.r')

# Fast function to get neighbor distances
# Refer to http://gis.stackexchange.com/questions/132384/distance-to-nearest-point-for-every-point-same-spatialpointsdataframe-in-r
getNeighbors <- function(dat, radius) {
  library(spdep)
  idlist <- dnearneigh(coordinates(dat), 0, radius)
  distlist <- nbdists(idlist, coordinates(dat))
  dflist <- list()
  pb <- txtProgressBar(0, length(idlist), style = 3)
  for (i in 1:length(idlist)) {
    if (any(distlist[[i]] <= radius)) {
      dflist[[i]] <- data.frame(idx = idlist[[i]], dist = distlist[[i]][distlist[[i]] <= radius])
    }
    else {
      dflist[[i]] <- NA
    }
	setTxtProgressBar(pb, i)
  }
  close(pb)
  dflist
}


radii <- c(1000,5000,7500,10000,20000,50000,100000)
fianhb_r <- getNeighbors(fiaalbers, radius = 5e5)

#save(fianhb_r, fiadist, trydist, pnwphylo, problemspp, fiataxa, fiaplotmat, fiaalbers, fiasums_plot, sppids, file = '/mnt/research/nasabio/data/fia/fiaworkspace.r')

save(fianhb_r, fiadist, trydist, pnwphylo, problemspp, fiataxa, fiaspatial, fiaalbers, fiacoords, fiasums_plot, sppids, fiaplotmat, file = '/mnt/research/nasabio/data/fia/fiaworkspace2.r')