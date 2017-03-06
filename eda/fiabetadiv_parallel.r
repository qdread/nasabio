# FIA taxonomic and phylogenetic beta diversity calculation for the cluster.
# Modified 6 March 2017: correct the bug in species names. Also add functional beta-diversity to this.

fp <- '/mnt/research/nasabio'
fiapnw <- read.csv(file.path(fp, 'data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)

#phylogenetic distance matrix
library(ape)
load(file.path(fp, 'data/fia/phytophylo_fia.r'))
fiadist <- cophenetic(fiaphytophylo)

# Get FIA species code to scientific name lookup table
fiataxa <- read.csv(file.path(fp, 'data/fia/fia_taxon_lookuptable.csv'), stringsAsFactors = FALSE)
pnw_codes <- unique(fiapnw$SPCD)
#all(pnw_codes %in% fiataxa$FIA.Code) #yes

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
#sppnames <- pnw_species$sciname[match(sppids, pnw_codes)]
dimnames(fiaplotmat)[[2]] <- fiataxa$sciname[idx]

# Get rid of the unknown species.
fiaplotmat <- fiaplotmat[, dimnames(fiaplotmat)[[2]] %in% fiaphytophylo$tip.label]

# Generate distance matrix for functional beta-diversity
source('~/code/fia/trydistmat.r')



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


radii <- c(1000,2000,3000,4000,5000,10000)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
r <- radii[task]


# Initialize data structures for observed metrics
fia_shannonbetadiv <- rep(NA, nrow = nrow(fiaalbers))
fia_meanpairwisedissim <- rep(NA, nrow = nrow(fiaalbers))
fia_meanpairwisedissim_pa <- rep(NA, nrow = nrow(fiaalbers))
fia_phypairwise <- rep(NA, nrow = nrow(fiaalbers))
fia_phypairwise_pa <- rep(NA, nrow = nrow(fiaalbers))
fia_phynt <- rep(NA, nrow = nrow(fiaalbers))
fia_phynt_pa <- rep(NA, nrow = nrow(fiaalbers))
fia_phypairwise_z <- rep(NA, nrow = nrow(fiaalbers))
fia_phypairwise_pa_z <- rep(NA, nrow = nrow(fiaalbers))
fia_phynt_z <- rep(NA, nrow = nrow(fiaalbers))
fia_phynt_pa_z <- rep(NA, nrow = nrow(fiaalbers))
fia_funcpairwise <- rep(NA, nrow = nrow(fiaalbers))
fia_funcpairwise_pa <- rep(NA, nrow = nrow(fiaalbers))
fia_funcnt <- rep(NA, nrow = nrow(fiaalbers))
fia_funcnt_pa <- rep(NA, nrow = nrow(fiaalbers))
fia_funcpairwise_z <- rep(NA, nrow = nrow(fiaalbers))
fia_funcpairwise_pa_z <- rep(NA, nrow = nrow(fiaalbers))
fia_funcnt_z <- rep(NA, nrow = nrow(fiaalbers))
fia_funcnt_pa_z <- rep(NA, nrow = nrow(fiaalbers))


fia_nneighb <- rep(NA, nrow = nrow(fiaalbers))


pb2 <- txtProgressBar(0, nrow(fiaalbers), style = 3)

# Number of simulations for null model
nnull <- 999

library(vegan)
library(vegetarian)

source('~/code/fia/fixpicante.r')
trydist <- as.matrix(trydist)

fianhb_r <- getNeighbors(fiaalbers, radius = r)
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
        mat_p_noproblem <- mat_p[, !dimnames(mat_p)[[2]] %in% problemspp, drop = FALSE]
        
        # Calculate beta-diversity for that matrix.
        
        fia_shannonbetadiv[p] <- d(abundances = mat_p, lev = 'beta', wts = FALSE, q = 1)
        fia_meanpairwisedissim[p] <- mean(vegdist(x = mat_p, binary = FALSE, method = 'jaccard'))
        fia_meanpairwisedissim_pa[p] <- mean(vegdist(x = mat_p, binary = TRUE, method = 'jaccard'))
        
        # Must catch errors in comdist for when the number of columns is zero
        if (ncol(mat_p) > 1) {
          fia_phypairwise[p] <- mean(comdist(comm = mat_p, dis = fiadist, abundance.weighted = TRUE))
          fia_phypairwise_pa[p] <- mean(comdist(comm = mat_p, dis = fiadist, abundance.weighted = FALSE))
          fia_phynt[p] <- mean(comdistnt(comm = mat_p, dis = fiadist, abundance.weighted = TRUE))
          fia_phynt_pa[p] <- mean(comdistnt(comm = mat_p, dis = fiadist, abundance.weighted = FALSE))
		  fia_funcpairwise[p] <- mean(comdist(comm = mat_p_noproblem, dis = trydist, abundance.weighted = TRUE))
          fia_funcpairwise_pa[p] <- mean(comdist(comm = mat_p_noproblem, dis = trydist, abundance.weighted = FALSE))
          fia_funcnt[p] <- mean(comdistnt(comm = mat_p_noproblem, dis = trydist, abundance.weighted = TRUE))
          fia_funcnt_pa[p] <- mean(comdistnt(comm = mat_p_noproblem, dis = trydist, abundance.weighted = FALSE))
		  
          # Null models by scrambling distance matrix
          phypairwise_null <- phypairwise_pa_null <- phynt_null <- phynt_pa_null <- rep(NA, nnull)
		  funcpairwise_null <- funcpairwise_pa_null <- funcnt_null <- funcnt_pa_null <- rep(NA, nnull)
          
          for (sim in 1:nnull) {
            nullidx <- sample(1:nrow(fiadist))
            fiadistnull <- fiadist
            dimnames(fiadistnull)[[1]] <- dimnames(fiadistnull)[[1]][nullidx]
            dimnames(fiadistnull)[[2]] <- dimnames(fiadistnull)[[2]][nullidx]
            
			trynullidx <- sample(1:nrow(trydist))
            trydistnull <- trydist
            dimnames(trydistnull)[[1]] <- dimnames(trydistnull)[[1]][trynullidx]
            dimnames(trydistnull)[[2]] <- dimnames(trydistnull)[[2]][trynullidx]
			
            phypairwise_null[sim] <- mean(comdist(comm = mat_p, dis = fiadistnull, abundance.weighted = TRUE))
            phypairwise_pa_null[sim] <- mean(comdist(comm = mat_p, dis = fiadistnull, abundance.weighted = FALSE))
            phynt_null[sim] <- mean(comdistnt(comm = mat_p, dis = fiadistnull, abundance.weighted = TRUE))
            phynt_pa_null[sim] <- mean(comdistnt(comm = mat_p, dis = fiadistnull, abundance.weighted = FALSE))
			
			funcpairwise_null[sim] <- mean(comdist(comm = mat_p, dis = trydistnull, abundance.weighted = TRUE))
            funcpairwise_pa_null[sim] <- mean(comdist(comm = mat_p, dis = trydistnull, abundance.weighted = FALSE))
            funcnt_null[sim] <- mean(comdistnt(comm = mat_p, dis = trydistnull, abundance.weighted = TRUE))
            funcnt_pa_null[sim] <- mean(comdistnt(comm = mat_p, dis = trydistnull, abundance.weighted = FALSE))
          }
          
          fia_phypairwise_z[p] <- (fia_phypairwise[p] - mean(phypairwise_null, na.rm=T))/sd(phypairwise_null, na.rm=T)
          fia_phypairwise_pa_z[p] <- (fia_phypairwise_pa[p] - mean(phypairwise_pa_null, na.rm=T))/sd(phypairwise_pa_null, na.rm=T)
          fia_phynt_z[p] <- (fia_phynt[p] - mean(phynt_null, na.rm=T))/sd(phynt_null, na.rm=T)
          fia_phynt_pa_z[p] <- (fia_phynt_pa[p] - mean(phynt_pa_null, na.rm=T))/sd(phynt_pa_null, na.rm=T)
		  
		  fia_funcpairwise_z[p] <- (fia_funcpairwise[p] - mean(funcpairwise_null, na.rm=T))/sd(funcpairwise_null, na.rm=T)
          fia_funcpairwise_pa_z[p] <- (fia_funcpairwise_pa[p] - mean(funcpairwise_pa_null, na.rm=T))/sd(funcpairwise_pa_null, na.rm=T)
          fia_funcnt_z[p] <- (fia_funcnt[p] - mean(funcnt_null, na.rm=T))/sd(funcnt_null, na.rm=T)
          fia_funcnt_pa_z[p] <- (fia_funcnt_pa[p] - mean(funcnt_pa_null, na.rm=T))/sd(funcnt_pa_null, na.rm=T)
        }
        else {
          fia_phypairwise[p] <- 0
          fia_phypairwise_pa[p] <- 0
          fia_phynt[p] <- 0
          fia_phynt_pa[p] <- 0
          fia_phypairwise_z[p] <- 0
          fia_phypairwise_pa_z[p] <- 0
          fia_phynt_z[p] <- 0
          fia_phynt_pa_z[p] <- 0
		  fia_funcpairwise[p] <- 0
          fia_funcpairwise_pa[p] <- 0
          fia_funcnt[p] <- 0
          fia_funcnt_pa[p] <- 0
          fia_funcpairwise_z[p] <- 0
          fia_funcpairwise_pa_z[p] <- 0
          fia_funcnt_z[p] <- 0
          fia_funcnt_pa_z[p] <- 0
        }
        fia_nneighb[p] <- nrow(mat_p) - 1
        
        
      }
    }
  }
  setTxtProgressBar(pb2, p)
}

close(pb2)

# Compile all of these values into a single data frame and save.
fia_betadiv <- data.frame(nneighb = fia_nneighb,
						  beta_td_shannon = fia_shannonbetadiv, 
						  beta_td_pairwise = fia_meanpairwisedissim, 
						  beta_pd_pairwise_z = fia_phypairwise_z, 
						  beta_pd_nearesttaxon_z = fia_phynt_z,
						  beta_pd_pairwise_presence_z = fia_phypairwise_pa_z,
						  beta_pd_nearesttaxon_presence_z = fia_phynt_pa_z,
						  beta_fd_pairwise_z = fia_funcpairwise_z,
						  beta_fd_nearesttaxon_z = fia_funcnt_z,
						  beta_fd_pairwise_presence_z = fia_funcpairwise_pa_z,
						  beta_fd_nearesttaxon_presence_z = fia_funcnt_pa_z)

write.csv(fia_betadiv, file = file.path(fp, paste0('fia_allbetadiv',task,'.csv')), row.names = FALSE)						  
