# FIA taxonomic and phylogenetic beta diversity calculation for the cluster.
# Forked 6 April 2017: Just do taxonomic diversity for now so that it will run faster.
# Modified 6 March 2017: correct the bug in species names. Also add functional beta-diversity to this.

fp <- '/mnt/research/nasabio'
fiapnw <- read.csv(file.path(fp, 'data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)

#phylogenetic distance matrix
library(ape)
load(file.path(fp, 'data/fia/phytophylo_fia.r'))
#fiadist <- cophenetic(fiaphytophylo)

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
fiataxa$sciname[fiataxa$sciname == 'Populus_balsamiferassp.Trichocarpa'] <- 'Populus_trichocarpa'


fiaplotlist <- fiasums_plot %>% do(x = area_by_sp(., sppids))
fiaplotmat <- do.call('rbind', fiaplotlist$x)

# Get the species names from the lookup table that go with the numerical codes.
dimnames(fiaplotmat)[[2]] <- fiataxa$sciname[idx]

# Get rid of the unknown species.
fiaplotmat <- fiaplotmat[, dimnames(fiaplotmat)[[2]] %in% fiaphytophylo$tip.label]

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


radii <- c(1000,5000,7500,10000,20000,50000)


# Initialize data structures for observed metrics
fia_shannonbetadiv <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_meanpairwisedissim <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))
fia_meanpairwisedissim_pa <- matrix(NA, nrow = nrow(fiaalbers), ncol = length(radii))

fia_nneighb <- matrix(0, nrow = nrow(fiaalbers), ncol = length(radii))

library(vegan)
library(vegetarian)

for (r in 1:length(radii)) {

pb2 <- txtProgressBar(0, nrow(fiaalbers), style = 3)

fianhb_r <- getNeighbors(fiaalbers, radius = radii[r])
for (p in 1:nrow(fiaalbers)) {
  if (class(fianhb_r[[p]]) == 'data.frame') {
    # Subset out the data frame with the nearest neighbors
    plotcns <- fiaalbers[c(p, fianhb_r[[p]]$idx), ]$PLT_CN
    dat_p <- subset(fiasums_plot, PLT_CN %in% plotcns)
    # Convert into a site x species matrix
    sppids <- sort(unique(dat_p$SPCD))
    mat_p <- dat_p %>% group_by(PLT_CN) %>% do(x = area_by_sp(., sppids))
    mat_p <- do.call('rbind', mat_p$x)
    
    if(!is.null(mat_p)) {
      if(nrow(mat_p) > 1) {
        # Fix the species names to match the phylogeny, and get rid of the unknown species.
        sppnames <- fiataxa$sciname[match(sppids, fiataxa$FIA.Code)]
        dimnames(mat_p)[[1]] <- 1:nrow(mat_p)
        dimnames(mat_p)[[2]] <- sppnames
        mat_p <- mat_p[, dimnames(mat_p)[[2]] %in% fiaphytophylo$tip.label, drop = FALSE]
        
        # Calculate beta-diversity for that matrix.
        
        fia_shannonbetadiv[p, r] <- d(abundances = mat_p, lev = 'beta', wts = FALSE, q = 1)
        fia_meanpairwisedissim[p, r] <- mean(vegdist(x = mat_p, binary = FALSE, method = 'jaccard'))
        fia_meanpairwisedissim_pa[p, r] <- mean(vegdist(x = mat_p, binary = TRUE, method = 'jaccard'))
        
        
        fia_nneighb[p, r] <- nrow(mat_p) - 1
        
        
      }
    }
  }
  setTxtProgressBar(pb2, p)
}

close(pb2)

}

# Compile all of these values into a single data frame and save.
m2v <- function(m) m[1:prod(dim(m))]

fia_betadiv <- cbind(ungroup(fiacoords),
				data.frame(radius = rep(radii, each=nrow(fiaalbers)),
					      nneighb = m2v(fia_nneighb),
						  beta_shannon = m2v(fia_shannonbetadiv), 
						  beta_pairwise_abundance = m2v(fia_meanpairwisedissim),
						  beta_pairwise_presence = m2v(fia_meanpairwisedissim_pa)
						  ))

write.csv(fia_betadiv, file = file.path(fp, 'data/fia/fia_betatd.csv'), row.names = FALSE)						  
