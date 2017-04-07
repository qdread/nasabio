# Create distance matrix from try data. Use for beta-diversity calculations.
# May need to tweak the community matrix to get rid of some of the species with no trait data.

library(FD)

try_wide <- read.csv(file.path(fp, 'data/fia/try_fia.csv'), stringsAsFactors = FALSE)

# Edit 07 Apr: fix to get the populus correct.
try_wide$AccSpeciesName[try_wide$AccSpeciesName == 'Populus trichocarpa'] <- 'Populus balsamifera trichocarpa'

tryspp <- gsub(' ', '_', try_wide$AccSpeciesName)
trymatch <- dimnames(fiaplotmat)[[2]] %in% tryspp
dimnames(fiaplotmat)[[2]][!trymatch]
# Aesculus sp and fraxinus sp don't match. Aesculus sp. can have the traits of Aesculus californica, Fraxinus sp can have average of the two Fraxinus in the try dataset.

try_use <- try_wide[,c(3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 20, 21, 22, 24, 25)]

aesculus <- try_use[tryspp == 'Aesculus_californica', ]
fraxinus <- try_use[grep('Fraxinus',tryspp), ]
fraxinus <- apply(fraxinus, 2, mean, na.rm=T)
fraxinus <- t(fraxinus)

tryspp <- c(tryspp, 'Aesculus_sp', 'Fraxinus_sp')
try_use <- rbind(try_use, aesculus, fraxinus)
dimnames(try_use)[[1]] <- tryspp

# Pare down try spp list to the ones in FIA, and get rid of species with too few traits.
try_use <- try_use[dimnames(fiaplotmat)[[2]], ]
ntraits <- apply(try_use, 1, function(x) sum(!is.na(x)))

# Check which plots have greater than 25% species that we don't have traits for. Throw them out later.
problemspp <- names(ntraits)[ntraits < 3]
percent_problem <- apply(fiaplotmat, 1, function(x) sum(x[names(x) %in% problemspp])/sum(x))
# Only a few hundred plots have greater than 25% so we can just throw them out.

fiaplotmat_noproblem <- fiaplotmat[, !dimnames(fiaplotmat)[[2]] %in% problemspp]
try_noproblem <- try_use[!dimnames(try_use)[[1]] %in% problemspp, ]

# Code also breaks if a plot has zero trees in it.
zerorows <- apply(fiaplotmat_noproblem, 1, sum) == 0

trydist <- gowdis(try_noproblem)