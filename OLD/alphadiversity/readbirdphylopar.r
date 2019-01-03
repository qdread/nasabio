# Read bird phylogenies and match to the real species names that we have. Also add species identified only to genus to the phylogeny.
# Author: QDR
# Project: Aquaxterra
# Created: 04 Nov 2016
# Last modified: 22 Feb 2017

# Modified 22 Feb: Non-abundance-weighted!
# Modified 08 Dec: Extensive debugging
# Modified 06 Dec: Use the updated BBS.
# Modified 05 Dec: instead of averaging the distance matrix, just make the whole thing repeat 10x with a single tree each time.
# Modified 03 Dec: average ten trees instead of using just one (follow Jarzyna); also get rid of nocturnal species.
# Modified 15 Nov: Get rid of years older than 1997 (when modern bbs data was established)
# Modified 10 Nov: Add PD calculations for one test tree
# Modified 08 Nov: Resolve AOUs at random for each tree to create a new tree.

# Extract task ID from the system.
task <- as.numeric(Sys.getenv('PBS_ARRAYID')) # Tasks 1-100 (10 trees x 10 matrix slices)

# Edit 12 Dec: redo the trees that didn't work the first time.
# Edited again on 13 Dec.
# Have to redo all the trees and slices for presence and residence.
tree_list <- rep(1:10, times=10)
slice_list <- rep(1:10, each=10)
tree_to_use <- tree_list[task]
slice_to_use <- slice_list[task]

#tree_to_use <- ceiling(task/10)
#slice_to_use <- rep(1:10, times = 10)[task]


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
# Edited 08 Dec: moved after, and get rid of the unused tips (should include all the NAs)

# ericsondist <- cophenetic(t1)
# Must sort the distance matrices so that they are all the same order.
# roworder <- dimnames(ericsondist[[1]])[[1]]
# ericsondist <- lapply(ericsondist, function(x) x[roworder, roworder])
# Mean branch lengths across the ten randomly sampled trees
# ericsondistmean <- apply(simplify2array(ericsondist), 1:2, mean)

# Quick correction to fix two birds that aren't in the phylogeny. Just get rid of the eastern yellow wagtail since it's probably only in Alaska anyway.
fixedbbsmat[, which(sppids == 5739)] <- fixedbbsmat[, which(sppids == 5738)] + fixedbbsmat[, which(sppids == 5739)]
fixedbbsmat[, which(sppids == 5738)] <- 0
fixedbbsmat[, which(sppids == 6960)] <- 0
#fixedbbsmat[, 'Artemisiospiza belli'] <- fixedbbsmat[, 'Artemisiospiza nevadensis'] + fixedbbsmat[, 'Artemisiospiza belli']
#fixedbbsmat[, 'Artemisiospiza nevadensis'] <- 0
#fixedbbsmat[, 'Motacilla tschutschensis'] <- 0

# Match the tip labels of ericson or hackett tree with the row names of the fixed bbs matrix.
ns <- colSums(fixedbbsmat)
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
dimnames(fixedbbsmat)[[2]] <- dimnames_matrix

fixedbbsmat_nonzero <- fixedbbsmat[, ns > 0]
fixedbbsmat_nonzero <- fixedbbsmat_nonzero[, !(dimnames(fixedbbsmat_nonzero)[[2]] %in% nocturnalbirds)]

t1 <- drop.tip(phy = t1, tip = which(!t1$tip.label %in% dimnames(fixedbbsmat_nonzero)[[2]]))
ericsondist <- cophenetic(t1)

# Actual calculation of pd, mpd, and mntd

library(picante)

# Split this into smaller jobs that can be run in parallel.


xx <- round(seq(0, nrow(fixedbbsmat), length.out=11))
xxmat <- cbind((xx+1)[-11], xx[-1])
rowstouse <- (xxmat[slice_to_use,1]:xxmat[slice_to_use,2])

x <- fixedbbsmat_nonzero[rowstouse,]

# Mean PD across ten trees
# allpds <- list()
# for (i in 1:10) {
	# allpds[[i]] <- pd(x, use_trees[[i]], include.root=TRUE)
# }
#pd_ericson <- apply(Reduce(cbind, lapply(allpds, function(x) x$PD)), 1, mean) # messed up way of getting average pd

pd_ericson <- pd(x, t1, include.root = TRUE)
mpd_ericson <- ses.mpd(x, ericsondist, null.model = 'independentswap', abundance.weighted = FALSE, runs = 999, iterations = 1000)
mntd_ericson <- ses.mntd(x, ericsondist, null.model = 'independentswap', abundance.weighted = FALSE, runs = 999, iterations = 1000)

save(pd_ericson, mpd_ericson, mntd_ericson, file = paste0('/mnt/research/aquaxterra/DATA/raw_data/bird_traits/bird_phylogeny/pd_tree_presence_',tree_to_use,'_part',slice_to_use, '.r'))
