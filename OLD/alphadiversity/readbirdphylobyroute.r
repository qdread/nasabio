# Read bird phylogenies and match to the real species names that we have. Also add species identified only to genus to the phylogeny.
# Aggregated by route.
# Author: QDR
# Project: Aquaxterra
# Created: 04 Nov 2016
# Last modified: 22 Feb 2017

# Modified 22 Feb: Non-abundance-weighted
# Modified 08 Dec: Extensive debugging
# Modified 06 Dec: Use the updated BBS.
# Modified 05 Dec: instead of averaging the distance matrix, just make the whole thing repeat 10x with a single tree each time.
# 04 Dec: Copied this from readbirdphylopar.r
# Modified 03 Dec: average ten trees instead of using just one (follow Jarzyna); also get rid of nocturnal species.
# Modified 15 Nov: Get rid of years older than 1997 (when modern bbs data was established)
# Modified 10 Nov: Add PD calculations for one test tree
# Modified 08 Nov: Resolve AOUs at random for each tree to create a new tree.

# Extract task ID from the system.
task <- as.numeric(Sys.getenv('PBS_ARRAYID')) # Only 10 tasks (1 for each tree). No need to slice matrix up since there aren't too many rows.
tree_to_use <- task


library(ape)

# Read the randomly sampled 1000-tree subsets of the birdtree.org phylogenies (Ericson and Hackett versions)

erictree <- read.nexus('/mnt/research/aquaxterra/DATA/raw_data/bird_traits/bird_phylogeny/ericson1000.tre')
#hacktree <- read.nexus('/mnt/research/aquaxterra/DATA/raw_data/bird_traits/bird_phylogeny/hackett1000.tre')

# Load one of ten randomly selected trees (random numbers generated on 05 Dec)
treeids <- c(716, 566, 377, 568, 977, 141, 104, 194, 944, 67)
t1 <- erictree[[treeids[tree_to_use]]]

# Find which species on our species list are not given as tips of the tree.

tlabel1 <- t1$tip.label
tlabel1 <- gsub('_', ' ', tlabel1)

bbsspp <- read.csv('/mnt/research/aquaxterra/DATA/raw_data/bird_traits/specieslist.csv', stringsAsFactors = FALSE)
load('/mnt/research/aquaxterra/DATA/raw_data/BBS/bbsmatconsolidated2015.r')

phymatch <- bbsspp$Latin_Name_clean %in% tlabel1 | bbsspp$Latin_Name_synonym %in% tlabel1 | bbsspp$Latin_Name_synonym2 %in% tlabel1
phymatchidx <- rep(NA,length(tlabel1))

for (i in 1:length(tlabel1)) {
	phymatchidx[i] <- c(which(bbsspp$Latin_Name_clean == tlabel1[i] | bbsspp$Latin_Name_synonym == tlabel1[i] | bbsspp$Latin_Name_synonym2 == tlabel1[i]), NA)[1]
}

# We must add some non-identified species to the tree.
# The protocol should be: if it is unknown between 2 or 3 species, assign the unknown individual randomly to one of those species
# If it is unknown in a genus, assign randomly to any species in that genus
# If it is unknown in a family, assign randomly to any species in that family
# If it is a hybrid, assign randomly to one of the parent species
# If it is a subspecies, it gets the same ID as the undifferentiated species.

#library(phytools)
#spp_not_in_tree <- bbsspp[!phymatch, ]
# Manually edit these names then reload
#write.csv(spp_not_in_tree, file='DATA/raw_data/bird_traits/sppnotinphylo.csv', row.names=FALSE)

# function to get pd out of each row

# Exclude nocturnal birds
bt <- read.csv('/mnt/research/aquaxterra/DATA/raw_data/bird_traits/birdtraitmerged.csv', stringsAsFactors=FALSE)
nocturnalbirds <- bt$Latin_Name[bt$Nocturnal == 1]

#fixedbbsmat <- fixedbbsmat[bbsgrps$year >= 1997, ]

# The mean distance matrix must be calculated before matching up all the dimnames, because otherwise NA row names are made.
# Edit 08 Dec: Move this below with a call to drop.tip (remove all unused tips which includes all the NAs)

#ericsondist <- cophenetic(t1)
# Must sort the distance matrices so that they are all the same order.
# roworder <- dimnames(ericsondist[[1]])[[1]]
# ericsondist <- lapply(ericsondist, function(x) x[roworder, roworder])
# Mean branch lengths across the ten randomly sampled trees
# ericsondistmean <- apply(simplify2array(ericsondist), 1:2, mean)

# Quick correction to fix two birds that aren't in the phylogeny. Just get rid of the eastern yellow wagtail since it's probably only in Alaska anyway.
fixedbbsmat_byroute[, which(sppids == 5739)] <- fixedbbsmat_byroute[, which(sppids == 5738)] + fixedbbsmat_byroute[, which(sppids == 5739)]
fixedbbsmat_byroute[, which(sppids == 5738)] <- 0
fixedbbsmat_byroute[, which(sppids == 6960)] <- 0

# Match the tip labels of ericson or hackett tree with the row names of the fixed bbs matrix.
ns <- colSums(fixedbbsmat_byroute)
tlabelaou <- bbsspp$AOU[phymatchidx]
aoustoadd <- sppids[!sppids %in% tlabelaou & ns > 0] # AOUs that need to be added

# Set dimnames of fixedbbsmat and tip labels of erictree to be the same.
dimnames_tlabel <- rep(NA, length(tlabelaou))
for (i in 1:length(tlabelaou)) {
	names_i <- bbsspp[bbsspp$AOU == tlabelaou[i], c('Latin_Name')]
	if (length(names_i)>0) dimnames_tlabel[i] <- names_i[1]
}
dimnames_matrix <- rep(NA, length(sppids))
for (i in 1:length(sppids)) {
	names_i <- bbsspp[bbsspp$AOU == sppids[i], c('Latin_Name')]
	if (length(names_i)>0) dimnames_matrix[i] <- names_i[1]
}

t1$tip.label <- dimnames_tlabel 
dimnames(fixedbbsmat_byroute)[[2]] <- dimnames_matrix

fixedbbsmat_byroute_nonzero <- fixedbbsmat_byroute[, ns > 0]
fixedbbsmat_byroute_nonzero <- fixedbbsmat_byroute_nonzero[, !(dimnames(fixedbbsmat_byroute_nonzero)[[2]] %in% nocturnalbirds)]

t1 <- drop.tip(phy = t1, tip = which(!t1$tip.label %in% dimnames(fixedbbsmat_byroute_nonzero)[[2]]))
ericsondist <- cophenetic(t1)

# Actual calculation of pd, mpd, and mntd

library(picante)

x <- fixedbbsmat_byroute_nonzero

pd_ericson <- pd(x, t1, include.root = TRUE)
mpd_ericson <- ses.mpd(x, ericsondist, null.model = 'independentswap', abundance.weighted = FALSE, runs = 999, iterations = 1000)
mntd_ericson <- ses.mntd(x, ericsondist, null.model = 'independentswap', abundance.weighted = FALSE, runs = 999, iterations = 1000)

save(pd_ericson, mpd_ericson, mntd_ericson, file = paste0('/mnt/research/aquaxterra/DATA/raw_data/bird_traits/bird_phylogeny/pd_tree_presence_',tree_to_use,'_byroute.r'))
