# Beta-diversity calculation for breeding bird survey data
# Run in parallel on HPCC (array job with a different radius for each task)
# QDR, 07 Mar 2017
# Project: NASABioXGeo

# Note to self: once the hpcc is available again, make sure that all this code matches the most recent code on there.

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

bbsalbers <- bbsalbers[has_coords, ]
bbsgrps <- bbsgrps[has_coords, ]

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
fixedbbsmat <- fixedbbsmat[has_coords, !(dimnames(fixedbbsmat)[[2]] %in% nocturnalbirds)]

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
      dflist[[i]] <- data.frame(idx = idlist[[i]], dist = distlist[[i]][distlist[[i]] <= radius])
    }
    else {
      dflist[[i]] <- NA
    }
  }
  dflist
}

radii <- c(1000,2000,3000,4000,5000,10000,15000,20000) # Try these for now.
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
r <- radii[task]

library(dplyr)

bbscov <- cbind(bbsgrps, bbsalbers)
bbsnhb_list <- bbscov %>% group_by(year, rteNo) %>% do(l = getNeighbors(., radius = r)) # Takes approx 8min on hpcc
# Flatten this into one list
bbsnhb_r <- do.call('c', bbsnhb_list$l)

# 7. Run all the presence-only metrics with appropriate null models using d, vegdist, and comdist.

# Initialize data structures for observed metrics (presence only)
bbs_meanpairwisedissim_pa <- rep(NA, nrow = nrow(bbsalbers))
bbs_phypairwise_pa <- rep(NA, nrow = nrow(bbsalbers))
bbs_phynt_pa <- rep(NA, nrow = nrow(bbsalbers))
bbs_phypairwise_pa_z <- rep(NA, nrow = nrow(bbsalbers))
bbs_phynt_pa_z <- rep(NA, nrow = nrow(bbsalbers))
bbs_funcpairwise_pa <- rep(NA, nrow = nrow(bbsalbers))
bbs_funcnt_pa <- rep(NA, nrow = nrow(bbsalbers))
bbs_funcpairwise_pa_z <- rep(NA, nrow = nrow(bbsalbers))
bbs_funcnt_pa_z <- rep(NA, nrow = nrow(bbsalbers))

bbs_nneighb <- rep(NA, nrow = nrow(bbsalbers))

pb2 <- txtProgressBar(0, nrow(bbsalbers), style = 3)

# Number of simulations for null model
nnull <- 999

library(vegan)
library(vegetarian)

source('~/code/fia/fixpicante.r')

for (p in 1:nrow(bbsalbers)) {
  if (class(bbsnhb_r[[p]]) == 'data.frame') {
    # Subset out the data frame with the nearest neighbors
	##### EDIT THIS
	year_p <- bbscov$year[p]
	route_p <- bbscov$rteNo[p]
	
    rowidx <- c(p, 
    dat_p <- fixedbbsmat[blablaindex, ]
    # Convert into a site x species matrix
	##### EDIT THIS
    sppids <- sort(unique(dat_p$SPCD))
    mat_p <- dat_p %>% group_by(PLT_CN) %>% do(x = area_by_sp(., sppids))
    mat_p <- do.call('rbind', mat_p$x)
    
    if(!is.null(mat_p)) {
      if(nrow(mat_p) > 1) {
		
		# EDIT THE BELOW 5 LINES!!!!!1
        # Fix the species names to match the phylogeny, and get rid of the unknown species.
        sppnames <- fiataxa$sciname[match(sppids, fiataxa$FIA.Code)]
        dimnames(mat_p)[[1]] <- 1:nrow(mat_p)
        dimnames(mat_p)[[2]] <- sppnames
        mat_p <- mat_p[, dimnames(mat_p)[[2]] %in% fiaphytophylo$tip.label, drop = FALSE]
        mat_p_noproblem <- mat_p[, !dimnames(mat_p)[[2]] %in% problemspp, drop = FALSE]
        
        # Calculate beta-diversity for that matrix.
        
        bbs_meanpairwisedissim_pa[p] <- mean(vegdist(x = mat_p, binary = TRUE, method = 'jaccard'))
        
        # Must catch errors in comdist for when the number of columns is zero
        if (ncol(mat_p) > 1) {
          bbs_phypairwise_pa[p] <- mean(comdist(comm = mat_p, dis = ericdist, abundance.weighted = FALSE))
          bbs_phynt_pa[p] <- mean(comdistnt(comm = mat_p, dis = ericdist, abundance.weighted = FALSE))
          bbs_funcpairwise_pa[p] <- mean(comdist(comm = mat_p, dis = birdtraitdist, abundance.weighted = FALSE))
          bbs_funcnt_pa[p] <- mean(comdistnt(comm = mat_p, dis = birdtraitdist, abundance.weighted = FALSE))
		  
          # Null models by scrambling distance matrix
          phypairwise_pa_null <- phynt_pa_null <- rep(NA, nnull)
		  funcpairwise_pa_null <- funcnt_pa_null <- rep(NA, nnull)
          
          for (sim in 1:nnull) {
            ericnullidx <- sample(1:nrow(ericdist))
            ericdistnull <- ericdist
            dimnames(ericdistnull)[[1]] <- dimnames(ericdistnull)[[1]][nullidx]
            dimnames(ericdistnull)[[2]] <- dimnames(ericdistnull)[[2]][nullidx]
            
			birdtraitnullidx <- sample(1:nrow(birdtraitdist))
            birdtraitdistnull <- birdtraitdist
            dimnames(birdtraitdistnull)[[1]] <- dimnames(birdtraitdistnull)[[1]][birdtraitnullidx]
            dimnames(birdtraitdistnull)[[2]] <- dimnames(birdtraitdistnull)[[2]][birdtraitnullidx]
			
            phypairwise_pa_null[sim] <- mean(comdist(comm = mat_p, dis = ericdistnull, abundance.weighted = FALSE))
            phynt_pa_null[sim] <- mean(comdistnt(comm = mat_p, dis = ericdistnull, abundance.weighted = FALSE))
			
            funcpairwise_pa_null[sim] <- mean(comdist(comm = mat_p, dis = birdtraitdistnull, abundance.weighted = FALSE))
            funcnt_pa_null[sim] <- mean(comdistnt(comm = mat_p, dis = birdtraitdistnull, abundance.weighted = FALSE))
          }
          
          bbs_phypairwise_pa_z[p] <- (bbs_phypairwise_pa[p] - mean(phypairwise_pa_null, na.rm=T))/sd(phypairwise_pa_null, na.rm=T)
          bbs_phynt_pa_z[p] <- (bbs_phynt_pa[p] - mean(phynt_pa_null, na.rm=T))/sd(phynt_pa_null, na.rm=T)
		  
          bbs_funcpairwise_pa_z[p] <- (bbs_funcpairwise_pa[p] - mean(funcpairwise_pa_null, na.rm=T))/sd(funcpairwise_pa_null, na.rm=T)
          bbs_funcnt_pa_z[p] <- (bbs_funcnt_pa[p] - mean(funcnt_pa_null, na.rm=T))/sd(funcnt_pa_null, na.rm=T)
        }
        else {
          bbs_phypairwise_pa[p] <- 0
          bbs_phynt_pa[p] <- 0
          bbs_phypairwise_pa_z[p] <- 0
          bbs_phynt_pa_z[p] <- 0
          bbs_funcpairwise_pa[p] <- 0
          bbs_funcnt_pa[p] <- 0
          bbs_funcpairwise_pa_z[p] <- 0
          bbs_funcnt_pa_z[p] <- 0
        }
        bbs_nneighb[p] <- nrow(mat_p) - 1
        
        
      }
    }
  }
  setTxtProgressBar(pb2, p)
}

close(pb2)

# Compile all of these values into a single data frame and save.
bbs_betadiv <- data.frame(nneighb = bbs_nneighb,
						  beta_td_pairwise_presence = bbs_meanpairwisedissim_pa, 
						  beta_pd_pairwise_presence_z = bbs_phypairwise_pa_z,
						  beta_pd_nearesttaxon_presence_z = bbs_phynt_pa_z,
						  beta_fd_pairwise_presence_z = bbs_funcpairwise_pa_z,
						  beta_fd_nearesttaxon_presence_z = bbs_funcnt_pa_z)

write.csv(bbs_betadiv, file = paste0('/mnt/research/nasabio/data/bbs/bbs_allbetadiv',task,'.csv'), row.names = FALSE)						  