# FIA taxonomic and phylogenetic beta diversity calculation for the cluster.

fp <- '/mnt/research/nasabio'
fiapnw <- read.csv(file.path(fp, 'data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)

library(ape)
load(file.path(fp, 'data/fia/phytophylo_fia.r'))
fiadist <- cophenetic(fiaphytophylo)

fiataxa <- read.csv(file.path(fp, 'data/fia/fia_taxon_lookuptable.csv'), stringsAsFactors = FALSE)
pnw_codes <- unique(fiapnw$SPCD)
all(pnw_codes %in% fiataxa$FIA.Code) #yes

# Scientific names of FIA species in PNW plots.
pnw_species <- fiataxa[match(pnw_codes, fiataxa$FIA.Code), c('Genus','Species')]
pnw_scinames <- paste(pnw_species$Genus, pnw_species$Species, sep = '_')
pnw_scinames <- gsub(' ', '', pnw_scinames) #Remove extraneous spaces

# Correction for subspecies that are not in the tree.
pnw_species$sciname <- pnw_scinames
pnw_species$sciname[pnw_species$sciname == 'Chrysolepis_chrysophyllavar.chrysophylla'] <- 'Chrysolepis_chrysophylla'
pnw_species$sciname[pnw_species$sciname == 'Populus_balsamiferassp.Trichocarpa'] <- 'Populus_trichocarpa'

pnw_species <- subset(pnw_species, !grepl('Tree', sciname)) # unidentified species
#pnw_species$sciname <- gsub('_', ' ', pnw_species$sciname)
pnw_species <- pnw_species[order(pnw_species$sciname), ]

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

fia_nneighb <- rep(NA, nrow = nrow(fiaalbers))


pb2 <- txtProgressBar(0, nrow(fiaalbers), style = 3)

# Number of simulations for null model
nnull <- 999

library(picante)
library(vegan)
library(vegetarian)

# Fix match.comm.dist to suppress output
match.comm.dist.fixed <- function (comm, dis) 
{
  res <- list()
  commtaxa <- colnames(comm)
  if (is.null(commtaxa)) {
    stop("Community data set lacks taxa (column) names, these are required to match distance matrix and community data")
  }
  disclass <- class(dis)
  dis <- as.matrix(dis)
  distaxa <- rownames(dis)
  if (is.null(distaxa)) {
    warning("Distance matrix lacks taxa names, these are required to match community and distance matrix. Data are returned unsorted. Assuming that distance matrix and community data taxa columns are in the same order!")
    if (disclass == "dist") {
      return(list(comm = comm, dist = as.dist(dis)))
    }
    else {
      return(list(comm = comm, dist = dis))
    }
  }
  if (!all(distaxa %in% commtaxa)) {
    #print("Dropping taxa from the distance matrix because they are not present in the community data:")
    #print(setdiff(distaxa, commtaxa))
    dis <- dis[intersect(distaxa, commtaxa), intersect(distaxa, 
                                                       commtaxa)]
    distaxa <- rownames(dis)
  }
  if (any(!(commtaxa %in% distaxa))) {
    #    print("Dropping taxa from the community because they are not present in the distance matrix:")
    #   print(setdiff(commtaxa, distaxa))
    res$comm <- comm[, intersect(commtaxa, distaxa)]
  }
  else {
    res$comm <- comm
  }
  if (disclass == "dist") {
    res$dist <- as.dist(dis[colnames(comm), colnames(comm)])
  }
  else {
    res$dist <- dis[colnames(comm), colnames(comm)]
  }
  return(res)
}
assignInNamespace('match.comm.dist', match.comm.dist.fixed, 'picante')

fianhb_r <- getNeighbors(fiaalbers, radius = r)
for (p in 1:nrow(fiaalbers)) {
  if (class(fianhb_r[[p]]) == 'data.frame') {
    # Subset out the data frame with the nearest neighbors
    plotcns <- fiaalbers[c(p, fianhb_r[[p]]$idx), ]$PLT_CN
    dat_p <- subset(fiasums_plot, PLT_CN %in% plotcns)
    # Convert into a site x species matrix
    sppids <- sort(unique(dat_p$SPCD))
    mat_p <- dat_p %>% group_by(PLT_CN) %>% do(x = sapply(sppids, function(z) sum(.$basalarea[.$SPCD == z])))
    mat_p <- do.call('rbind', mat_p$x)
    
    if(!is.null(mat_p)) {
      if(nrow(mat_p) > 1) {
        # Fix the species names to match the phylogeny, and get rid of the unknown species.
        sppnames <- pnw_species$sciname[match(sppids, pnw_codes)]
        dimnames(mat_p)[[1]] <- 1:nrow(mat_p)
        dimnames(mat_p)[[2]] <- sppnames
        mat_p <- mat_p[, dimnames(mat_p)[[2]] %in% fiaphytophylo$tip.label, drop = FALSE]
        
        
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
          # Null models by scrambling distance matrix
          phypairwise_null <- phypairwise_pa_null <- phynt_null <- phynt_pa_null <- rep(NA, nnull)
          
          for (sim in 1:nnull) {
            nullidx <- sample(1:nrow(fiadist))
            fiadistnull <- fiadist
            dimnames(fiadistnull)[[1]] <- dimnames(fiadistnull)[[1]][nullidx]
            dimnames(fiadistnull)[[2]] <- dimnames(fiadistnull)[[2]][nullidx]
            
            phypairwise_null[sim] <- mean(comdist(comm = mat_p, dis = fiadistnull, abundance.weighted = TRUE))
            phypairwise_pa_null[sim] <- mean(comdist(comm = mat_p, dis = fiadistnull, abundance.weighted = FALSE))
            phynt_null[sim] <- mean(comdistnt(comm = mat_p, dis = fiadistnull, abundance.weighted = TRUE))
            phynt_pa_null[sim] <- mean(comdistnt(comm = mat_p, dis = fiadistnull, abundance.weighted = FALSE))
          }
          
          fia_phypairwise_z[p] <- (fia_phypairwise[p] - mean(phypairwise_null, na.rm=T))/sd(phypairwise_null, na.rm=T)
          fia_phypairwise_pa_z[p] <- (fia_phypairwise_pa[p] - mean(phypairwise_pa_null, na.rm=T))/sd(phypairwise_pa_null, na.rm=T)
          fia_phynt_z[p] <- (fia_phynt[p] - mean(phynt_null, na.rm=T))/sd(phynt_null, na.rm=T)
          fia_phynt_pa_z[p] <- (fia_phynt_pa[p] - mean(phynt_pa_null, na.rm=T))/sd(phynt_pa_null, na.rm=T)
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
        }
        fia_nneighb[p] <- nrow(mat_p) - 1
        
        
      }
    }
  }
  setTxtProgressBar(pb2, p)
}

close(pb2)

save(list = grep('fia_', ls(), value=TRUE), file = file.path(fp, paste0('fia_taxonomicbetadiv',task,'.RData')))
