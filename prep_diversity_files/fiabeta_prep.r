# FIA taxonomic and phylogenetic beta diversity calculation for the cluster.
# Forked 17 April. Just preparatory work, actual diversity calculated in a different script.
# Modified 17 May: add rows to matrix that have zero species.
# Modified 22 Aug: add better imputed traits
# Modified 20 Oct: new coordinates (workspace should now only contain pointers to the true coordinates)

fp <- '/mnt/research/nasabio'
fiapnw <- read.csv(file.path(fp, 'data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)

#phylogenetic distance matrix
library(ape)
load(file.path(fp, 'data/fia/pnwphylo_potter.r'))
fiadist <- cophenetic(pnwphylo)

# Get FIA species code to scientific name lookup table
fiataxa <- read.csv(file.path(fp, 'data/fia/fia_taxon_lookuptable.csv'), stringsAsFactors = FALSE)
pnw_codes <- unique(fiapnw$SPCD)

# Correct some subspecies IDs in fiataxa.
fiataxa$sciname <- paste(gsub(' ', '', fiataxa$Genus), gsub(' ', '', fiataxa$Species), sep = '_')
fiataxa$sciname[fiataxa$sciname == 'Chrysolepis_chrysophyllavar.chrysophylla'] <- 'Chrysolepis_chrysophylla'
fiataxa$sciname[fiataxa$sciname == 'Populus_balsamiferassp.Trichocarpa'] <- 'Populus_balsamifera_trichocarpa'

# Here, sum the Abies shastensis with Abies magnifica and get rid of Abies shastensis.
# Previously this was done too late in the pipeline. We don't want this currently invalid species cropping up later.
ashast <- fiataxa$FIA.Code[fiataxa$sciname == 'Abies_shastensis']
amagn <- fiataxa$FIA.Code[fiataxa$sciname == 'Abies_magnifica']

fiapnw$SPCD[fiapnw$SPCD == ashast] <- amagn			

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

fiacoords <- read.csv('/mnt/home/qdr/data/pnw.csv', stringsAsFactors = FALSE)
fiacoords <- fiacoords %>%
	filter(!is.na(ACTUAL_LAT)) %>%
	rename(PLT_CN = CN, lat = ACTUAL_LAT, lon = ACTUAL_LON)

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
traits_imputed <- read.csv(file.path(fp,'data/fia/traits_imputed_22aug.csv'), stringsAsFactors = FALSE, row.names = 1)
library(FD)

trydist <- gowdis(scale(traits_imputed))

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

fianhb_r <- getNeighbors(fiaalbers, radius = 5e5)

#save(fianhb_r, fiadist, trydist, traits_imputed, pnwphylo, fiataxa, fiaspatial, fiaalbers, fiacoords, fiasums_plot, sppids, fiaplotmat, file = '/mnt/research/nasabio/data/fia/fiaworkspace2.r')

# New saved workspace that does not actually contain the coordinates.
save(fianhb_r, fiadist, trydist, traits_imputed, pnwphylo, fiataxa, fiasums_plot, sppids, fiaplotmat, file = '/mnt/research/nasabio/data/fia/fiaworkspace_nospatial.r')
save(fiaspatial, fiaalbers, fiacoords, file = '/mnt/home/qdr/data/fiaworkspace_spatial.r')
write.csv(fiaalbers@data, file = '/mnt/research/nasabio/data/fia/fianocoords.csv', row.names = FALSE)