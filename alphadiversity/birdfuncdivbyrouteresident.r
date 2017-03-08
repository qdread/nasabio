# Bird functional diversity calculations: aggregated by route
# QDR/Aquaxterra/created 06Dec2016 (copied from birdfuncdivpar.r)

# Modified 22 Feb: Non-abundance-weighted
# Modified 07 Dec: Debug a lot.
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

ns <- colSums(fixedbbsmat_byroute)
fixedbbsmat_byroute_nonzero <- fixedbbsmat_byroute[, ns > 0 & !(sppids %in% nocturnalbirdAOUs) & !(sppids %in% migrantAOUs)]
sppids_nonzero <- sppids[ns > 0 & !(sppids %in% nocturnalbirdAOUs) & !(sppids %in% migrantAOUs)]

# Clean trait matrix and sort trait and sitexsp matrices so their dimensions all match.
dimnames(fixedbbsmat_byroute_nonzero)[[2]] <- sppids_nonzero # Already sorted by AOU
birdtraitclean <- birdtrait_diurnal[match(sppids_nonzero, birdtrait_diurnal$AOU), traitnames]
dimnames(birdtraitclean)[[1]] <- sppids_nonzero

# Make sure all columns are numerics, even the binary variables. (nocturnal is no longer included here)
birdtraitclean <- transform(birdtraitclean, PelagicSpecialist = as.numeric(PelagicSpecialist))

# Remove rows where no birds at all were found, flagging them for later. Probably none for entire routes.
rs <- apply(fixedbbsmat_byroute_nonzero, 1, sum)
fixedbbsmat_byroute_nonzerorows <- fixedbbsmat_byroute_nonzero[rs > 0, ]

# Run "naive" functional diversity, throwing all the traits in together.
library(FD)

# No need to parallelize if doing by route. there are only 50k route/year combos.

A <- fixedbbsmat_byroute_nonzerorows

zerocols <- apply(A, 2, sum) == 0

X <- birdtraitclean[!zerocols, ]

fd_all <- dbFD(x = X, a = A[, !zerocols], w.abun = FALSE, corr = 'cailliez')
save(fd_all, zerocols, file = file.path(fp, 'birdfuncdivobjectbyrouteresidentpresence.r'))
