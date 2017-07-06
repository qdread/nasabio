# Beta-diversity calculation for breeding bird survey data
# Run in parallel on HPCC (array job with a different radius for each task)
# QDR, 07 Mar 2017
# Project: NASABioXGeo

# Modification 11 April 2017: taking way too long so create an array job by radius.
# Reforked on 10 April 2017: Only preparatory work. Save "chunks" to send small pieces to HPCC.
# Forked on 06 April 2017: do taxonomic diversity only, which makes it quicker. We can just do 1 task per radius now.
# Modification 08 Mar 2017: Make doubly parallel by giving a different radius to each task
# but then also splitting the loop into 10 chunks since it will take forever to run otherwise.

radii <- c(1000,2000,3000,4000,5000,10000,15000,20000)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

combos <- data.frame(slice = rep(1:30, times = length(radii)), radius = rep(radii, each = 30))

r <- combos$radius[task]
slice <- combos$slice[task]

# 1. Define community either as aggregated route or by stop

# Start with by-stop diversity.

# 2. Load community matrix.

bbsspp <- read.csv('/mnt/research/aquaxterra/DATA/raw_data/bird_traits/specieslist.csv', stringsAsFactors = FALSE)
load('/mnt/research/aquaxterra/DATA/raw_data/BBS/bbsmatconsolidated2016.r') # Load fixed bbsmat. This loads both byroute and bystop.
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

#### edit on 5 july: make the matrix including the rows that don't have coords under the old way.
# fixedbbsmat <- fixedbbsmat[rs != 0, !(dimnames(fixedbbsmat)[[2]] %in% nocturnalbirds) & ns != 0]
# bbsalbers <- bbsalbers[rs != 0, ]
# bbsgrps <- bbsgrps[rs != 0, ]
# bbscov <- cbind(bbsgrps, bbsalbers)
# save(bbscov, fixedbbsmat, file = '/mnt/research/nasabio/data/bbs/bbsworkspace2016.r')

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

# Calculate distances and identities of all neighbors within the maximum radius
bbscov <- cbind(bbsgrps, bbsalbers)

# For optimization purposes, convert stop to numeric, and convert covariates to a matrix.
makenum <- function(x) {
  idx <- gregexpr('[0-9]', x)
  as.numeric(substr(x, idx[[1]][1], idx[[1]][length(idx[[1]])]))
}
bbscovmat <- bbscov %>% mutate(Stop = sapply(Stop, makenum), rteNo = as.numeric(rteNo)) %>% as.matrix

bbsnhb_list <- as.data.frame(bbscovmat) %>% group_by(year, rteNo) %>% do(l = getNeighbors(., radius = max(radii))) # Takes approx 8min on hpcc
# Flatten this into one list
bbsnhb_r <- do.call('c', bbsnhb_list$l)

save(bbsnhb_r, bbscov, bbscovmat, fixedbbsmat, file = '/mnt/research/nasabio/data/bbs/bbsworkspace.r')

