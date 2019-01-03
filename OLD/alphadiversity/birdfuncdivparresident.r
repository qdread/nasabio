# Bird functional diversity calculations
# QDR/Aquaxterra/created 17Nov2016/modified 06Dec2016

# Modified 07 Dec: Debug a lot.
# Modified 06 Dec: Use the updated BBS.
# Modified 05 Dec: Remove nocturnal birds (and remove nocturnal from the trait pca)

fp <- '/mnt/research/aquaxterra/DATA/raw_data/bird_traits'
birdtrait <- read.csv(file.path(fp, 'birdtraitmerged.csv'), stringsAsFactors = FALSE)

birdtrait[birdtrait == -999] <- NA
nocturnalbirdAOUs <- birdtrait$AOU[birdtrait$Nocturnal == 1]
migrantAOUs <- birdtrait$AOU[birdtrait$migrant_status == 1]
birdtrait_diurnal <- subset(birdtrait, (Nocturnal != 1 | is.na(Nocturnal)) & migrant_status == 0)

# Select traits to use
# Ones that seem important and/or have a lot of records.
# Includes diet, foraging strategy, body size, and life history.
traitnames <- names(birdtrait)[c(15:24, 29:36, 46:50, 53, 55:59)]

# Use the consolidated matrix that was used for the phylogenetic diversity calculations.
load('/mnt/research/aquaxterra/DATA/raw_data/BBS/bbsmatconsolidated2015.r') # Load fixed bbsmat.

# Small tweak added on Dec 9: AOUs 4451 and 3089 may be problematic. Get rid of them.

ns <- colSums(fixedbbsmat)
fixedbbsmat_nonzero <- fixedbbsmat[, ns > 0 & !(sppids %in% nocturnalbirdAOUs) & !(sppids %in% migrantAOUs) & !(sppids %in% c('3089','4451'))]
sppids_nonzero <- sppids[ns > 0 & !(sppids %in% nocturnalbirdAOUs) & !(sppids %in% migrantAOUs) & !(sppids %in% c('3089','4451'))]

# Clean trait matrix and sort trait and sitexsp matrices so their dimensions all match.
dimnames(fixedbbsmat_nonzero)[[2]] <- sppids_nonzero # Already sorted by AOU
birdtraitclean <- birdtrait_diurnal[match(sppids_nonzero, birdtrait_diurnal$AOU), traitnames]
dimnames(birdtraitclean)[[1]] <- sppids_nonzero

# Make sure all columns are numerics, even the binary variables. (nocturnal is no longer included here)
birdtraitclean <- transform(birdtraitclean, PelagicSpecialist = as.numeric(PelagicSpecialist))

# Remove rows where no birds at all were found, flagging them for later.
rs <- apply(fixedbbsmat_nonzero, 1, sum)
fixedbbsmat_nonzerorows <- fixedbbsmat_nonzero[rs > 0, ]

# Run "naive" functional diversity, throwing all the traits in together.
library(FD)

# Gower distance matrix.
# Mark pelagic specialist trait as a binary trait.
# The "Diet.*" and "ForStrat.*" traits are all percentages.
# The remaining traits are continuous traits with different units.
# Argument w can be added to weight traits differently, but currently not weighted.

# birdtrait_gowdist <- gowdis(x = birdtraitclean)

# Attempt to construct distance matrix without NAs by removing the three species with very high NA's (they are also taxonomically weird)
# This problem was fixed in another script.
#badIDs <- c('3151', '3410', '3460')

# As of 07 Dec, the bad AOUs are: 5738, 6960 (also fixed in another script)

#birdtraitgoodids <- birdtraitclean[!dimnames(birdtraitclean)[[1]] %in% badIDs, ]
# This works so we need to get rid of those bad IDs for the functional diversity. 

# Added 28 Nov: Since the code takes too long to run, split it into 10 bins and run it for each of the bins.

task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
rowidx <- round(seq(0,nrow(fixedbbsmat_nonzerorows),length.out=11))
rowidxmin <- rowidx[task]+1
rowidxmax <- rowidx[task+1]

A <- fixedbbsmat_nonzerorows[rowidxmin:rowidxmax, ]

zerocols <- apply(A, 2, sum) == 0

X <- birdtraitclean[!zerocols, ]

fd_all <- dbFD(x = X, a = A[, !zerocols], w.abun = FALSE, corr = 'cailliez')
save(fd_all, zerocols, file = file.path(fp, paste0('birdfuncdivobjectresidentpresence', task, '.r')))

# Fill back in the stops with zero birds, with zero for the species richness and NA for all functional diversity values.

# bfd_objects <- list()

# for (i in 1:10) {
	# load(file.path(fp, paste0('birdfuncdivobject', i, '.r')))
	# bfd_objects[[i]] <- with(fd_all, data.frame(nbsp, sing.sp, FRic, FEve, FDiv, FDis, RaoQ))
# }

# bfd_df <- do.call('rbind', bfd_objects)

# # Create df to hold the results, including the stops with zero birds.
# fd_all <- data.frame(bbsgrps, nbsp=0, sing.sp=0, FRic=NA, FEve=NA, FDiv=NA, FDis=NA, RaoQ=NA)
# fd_all[rs != 0, 4:10] <- bfd_df

# save(fd_all, file = file.path(fp, 'birdfuncdivobjectall.r'))
