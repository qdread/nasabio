# Beta-diversity calculation for breeding bird survey data
# Run in parallel on HPCC (array job with a different radius for each task)
# QDR, 07 Mar 2017
# Project: NASABioXGeo

# Forked on 06 April 2017: do taxonomic diversity only, which makes it quicker. We can just do 1 task per radius now.
# Modification 08 Mar 2017: Make doubly parallel by giving a different radius to each task
# but then also splitting the loop into 10 chunks since it will take forever to run otherwise.

radii <- c(1000,2000,3000,4000,5000,10000,15000,20000)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

r <- radii[task]

# 1. Define community either as aggregated route or by stop

# Start with by-stop diversity.

# 2. Load community matrix.

bbsspp <- read.csv('/mnt/research/aquaxterra/DATA/raw_data/bird_traits/specieslist.csv', stringsAsFactors = FALSE)
load('/mnt/research/aquaxterra/DATA/raw_data/BBS/bbsmatconsolidated2015.r') # Load fixed bbsmat. This loads both byroute and bystop.
# Quick correction to fix two birds that aren't in the phylogeny. Just get rid of the eastern yellow wagtail since it's probably only in Alaska anyway.
fixedbbsmat[, which(sppids == 5739)] <- fixedbbsmat[, which(sppids == 5738)] + fixedbbsmat[, which(sppids == 5739)]
fixedbbsmat[, which(sppids == 5738)] <- 0
fixedbbsmat[, which(sppids == 6960)] <- 0

# 3. Project lat-long to Albers. (Already done in another script, just load coordinates)

bbsalbers <- read.csv('/mnt/research/aquaxterra/DATA/raw_data/BBS/bbs_aea_coords.csv')

# 4. Create cophenetic distance matrix from bird phylogeny. Use average branch lengths from the ten trees.

library(ape)
erictree <- read.nexus('/mnt/research/aquaxterra/DATA/raw_data/bird_traits/bird_phylogeny/ericson1000.tre')
treeids <- c(716, 566, 377, 568, 977, 141, 104, 194, 944, 67)
# get distance matrix for each tree.
ericdistlist <- list()
for (i in treeids) ericdistlist[[length(ericdistlist) + 1]] <- cophenetic(erictree[[i]])

# sort distance matrices so that the species names in rows and columns are the same for each.
row_names <- dimnames(ericdistlist[[1]])[[1]]
column_names <- dimnames(ericdistlist[[1]])[[2]]
ericdistlist <- lapply(ericdistlist, function(x) x[row_names, column_names])

# take mean phylogenetic distance matrix
ericdist <- Reduce('+', ericdistlist)/length(ericdistlist)

# 5. Create trait distance matrix (Gower) from bird traits.

birdtrait <- read.csv('/mnt/research/aquaxterra/DATA/raw_data/bird_traits/birdtraitmerged.csv', stringsAsFactors = FALSE)

library(FD)

birdtrait[birdtrait == -999] <- NA
nocturnalbirdAOUs <- birdtrait$AOU[birdtrait$Nocturnal == 1]
birdtrait_diurnal <- subset(birdtrait, Nocturnal != 1 | is.na(Nocturnal))
traitnames <- names(birdtrait)[c(15:24, 29:36, 46:50, 53, 55:59)]

birdtraitdist <- as.matrix(gowdis(birdtrait_diurnal[, traitnames]))
dimnames(birdtraitdist)[[1]] <- dimnames(birdtraitdist)[[2]] <- birdtrait_diurnal$Latin_Name

# 5b. Perform any last data cleaning that is necessary.

# Remove stops that don't have coordinates from the dataset.
has_coords <- !is.na(bbsalbers$x_aea)

# Set dimnames of community matrix, trait distance matrix, and phylogenetic distance matrix to match.

dimnames_matrix <- rep(NA, length(sppids))
for (i in 1:length(sppids)) {
	names_i <- bbsspp[bbsspp$AOU == sppids[i], c('Latin_Name')]
	if (length(names_i)>0) dimnames_matrix[i] <- names_i[1]
}

dimnames(fixedbbsmat)[[2]] <- dimnames_matrix

tlabel1 <- erictree[[treeids[1]]]$tip.label
tlabel1 <- gsub('_', ' ', tlabel1)

phymatchidx <- rep(NA,length(tlabel1))

for (i in 1:length(tlabel1)) {
	phymatchidx[i] <- c(which(bbsspp$Latin_Name_clean == tlabel1[i] | bbsspp$Latin_Name_synonym == tlabel1[i] | bbsspp$Latin_Name_synonym2 == tlabel1[i]), NA)[1]
}
tlabelaou <- bbsspp$AOU[phymatchidx]

dimnames_tlabel <- rep(NA, length(tlabelaou))
for (i in 1:length(tlabelaou)) {
	names_i <- bbsspp[bbsspp$AOU == tlabelaou[i], c('Latin_Name')]
	if (length(names_i)>0) dimnames_tlabel[i] <- names_i[1]
}

dimnames(ericdist)[[1]] <- dimnames(ericdist)[[2]] <- dimnames_tlabel

# Remove nocturnal species, rows without coordinates, and rows and columns with zero sum.
ns <- colSums(fixedbbsmat)
rs <- rowSums(fixedbbsmat)
nocturnalbirds <- birdtrait$Latin_Name[birdtrait$Nocturnal == 1]
fixedbbsmat <- fixedbbsmat[has_coords & rs != 0, !(dimnames(fixedbbsmat)[[2]] %in% nocturnalbirds) & ns != 0]

bbsalbers <- bbsalbers[has_coords & rs != 0, ]
bbsgrps <- bbsgrps[has_coords & rs != 0, ]

# 6. For the given radius, get pairwise distances among stops and the list of all neighbors. Calculation can be sped up by limiting stop-level beta-diversity to stops that share a route.

# Fast function to get neighbor distances
# Refer to http://gis.stackexchange.com/questions/132384/distance-to-nearest-point-for-every-point-same-spatialpointsdataframe-in-r
getNeighbors <- function(dat, radius) {
  library(spdep)
  coords <- cbind(dat$x_aea, dat$y_aea)
  idlist <- dnearneigh(coords, 0, radius)
  distlist <- nbdists(idlist, coords)
  dflist <- list()
  for (i in 1:length(idlist)) {
    if (any(distlist[[i]] <= radius)) {
      dflist[[i]] <- data.frame(idx = idlist[[i]], Stop = dat$Stop[idlist[[i]]], dist = distlist[[i]][distlist[[i]] <= radius])
    }
    else {
      dflist[[i]] <- NA
    }
  }
  dflist
}

library(dplyr)

bbscov <- cbind(bbsgrps, bbsalbers)
bbsnhb_list <- bbscov %>% group_by(year, rteNo) %>% do(l = getNeighbors(., radius = r)) # Takes approx 8min on hpcc
# Flatten this into one list
bbsnhb_r <- do.call('c', bbsnhb_list$l)

# 7. Run all the presence-only metrics with appropriate null models using d, vegdist, and comdist.

# Initialize data structures for observed metrics (presence only)
#bbs_meanpairwisedissim_pa <- rep(NA, nrow(fixedbbsmat))
#bbs_nneighb <- rep(NA, nrow(fixedbbsmat))

library(vegan)
library(vegetarian)
library(foreach)
library(doParallel)

registerDoParallel(cores = 10)

# Done in parallel
bbs_list <- foreach(p = 1:nrow(fixedbbsmat)) %dopar% {
  if (class(bbsnhb_r[[p]]) == 'data.frame') {
    # Subset out the data frame with the nearest neighbors
	year_p <- bbscov$year[p]
	route_p <- bbscov$rteNo[p]
	stops_p <- c(as.character(bbscov$Stop[p]), as.character(bbsnhb_r[[p]]$Stop))
    dat_p <- fixedbbsmat[bbscov$year == year_p & bbscov$rteNo == route_p & bbscov$Stop %in% stops_p, ]
    # Get rid of empty columns
    mat_p <- dat_p[, apply(dat_p, 2, sum) > 0, drop = FALSE]
    
    if(!is.null(mat_p)) {
      if(nrow(mat_p) > 1) {
		
        dimnames(mat_p)[[1]] <- 1:nrow(mat_p)
        
        # Calculate beta-diversity for that matrix.
        #bbs_meanpairwisedissim_pa[p] <- mean(vegdist(x = mat_p, binary = TRUE, method = 'jaccard'))
        #bbs_nneighb[p] <- nrow(mat_p) - 1
        
		data.frame(nneighb = nrow(mat_p) - 1, beta_td_pairwise_presence = mean(vegdist(x = mat_p, binary = TRUE, method = 'jaccard')))
        
      }
	  else {
	    data.frame(nneighb = NA, beta_td_pairwise_presence = NA)
	  }
    }
	else {
	  data.frame(nneighb = NA, beta_td_pairwise_presence = NA)
	}
  }
  else {
    data.frame(nneighb = NA, beta_td_pairwise_presence = NA)
  }
}


# Compile all of these values into a single data frame and save.
bbs_betadiv <- do.call('rbind', bbs_list)

write.csv(bbs_betadiv, file = paste0('/mnt/research/nasabio/data/bbs/bbs_betatd',task,'.csv'), row.names = FALSE)						  