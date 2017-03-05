# Functional diversity of FIA
fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'

# First get list of species IDs for TRY.
tryspp <- read.csv(file.path(fp, 'tryspp.csv'), stringsAsFactors = FALSE)

trymatch <- pnw_species$sciname %in% tryspp$AccSpeciesName
pnw_species$sciname[!trymatch] # All match except for the "unknowns"

tryids <- tryspp$AccSpeciesID[match(pnw_species$sciname, tryspp$AccSpeciesName)]
write.table(t(na.omit(tryids)), sep=',')

# Load the try species list.
# Must disable quoting to read the entire file.
try_all <- read.delim(file.path(fp, 'nasa_try.txt'), stringsAsFactors = FALSE, quote = '')

nmeas <- table(try_all$AccSpeciesName, try_all$TraitName)
unittable <- table(try_all$TraitName, try_all$OrigUnitStr)

# Cut down try_all some more
try_all <- try_all[,c('DatasetID','AccSpeciesName','ObservationID','TraitName','DataName','OrigValueStr','UnitName','OrigUncertaintyStr','UncertaintyName')]

# Get rid of metadata rows for now
try_nometa <- subset(try_all, TraitName != "")

# Figure out whether individual traits have more than one unit of measurement.
measByUnitTable <- table(try_nometa$TraitName, try_nometa$UnitName)
measByUnitTable[apply(measByUnitTable>0, 1, sum) > 1, ]

# Plant longevity has some blank units and some in years
longevitynounit <- subset(try_nometa, grepl('Plant lifespan',TraitName) & UnitName=='')
#head(longevitynounit) # These are categorical longevity values (e.g. annual, perennial)
# Replace the names with two different values
try_nometa$TraitName[grepl('Plant lifespan',try_nometa$TraitName) & try_nometa$UnitName==''] <- 'Plant lifespan categorical'
try_nometa$TraitName[grepl('Plant lifespan',try_nometa$TraitName) & try_nometa$UnitName=='year'] <- 'Plant lifespan years'

# Seedbank longevity has some percentage, some dimensionless, and some in years
seedbanknounit <- subset(try_nometa, grepl('\\(seedbank\\) longevity',TraitName) & UnitName=='dimensionless') # Categorical (e.g. transient, persistent)
#seedbankpct <- subset(try_nometa, grepl('\\(seedbank\\) longevity',TraitName) & UnitName=='%') # Fraction of plots with persistent seeds

try_nometa$TraitName[grepl('\\(seedbank\\) longevity',try_nometa$TraitName) & try_nometa$UnitName=='dimensionless'] <- 'Seedbank longevity categorical'
#try_nometa$TraitName[grepl('\\(seedbank\\) longevity',try_nometa$TraitName) & try_nometa$UnitName=='%'] <- 'Seedbank persistence percent'

# Species by trait table for all the traits
spByTraitTable <- table(try_nometa$AccSpeciesName, try_nometa$TraitName)
apply(spByTraitTable>0, 2, sum) # How many species do we have for each trait?

library(reshape2)


mean_with_char <- function(x) {
  xnum <- as.numeric(x)
  if (any(!is.na(xnum))) as.character(mean(xnum, na.rm=TRUE)) else x[1]
}

# Function to change columns in a df that only consist of numeric strings to numerics.
numstring2num <- function(x) {
  xnum <- as.numeric(x)
  if (!any(is.na(xnum)) & !is.factor(x)) xnum else x
}

# Cast (convert from long form to wide so that each trait has its own column and can be edited)
try_byobsmean <- dcast(try_nometa[,c(1,2,3,4,5,6)], ObservationID+DatasetID+AccSpeciesName+DataName ~ TraitName, value.var='OrigValueStr', fun.aggregate = mean_with_char)

mean_and_category <- function(x) {
  xnum <- as.numeric(x$OrigValueStr)
  if (any(!is.na(xnum))) {
    data.frame(value = mean(xnum, na.rm=TRUE), category = NA)
  }
  else {
    data.frame(value = NA, category = x$OrigValueStr[1])
  }
}

library(dplyr)
try_byspmean <- try_nometa %>%
  group_by(AccSpeciesName, TraitName, DataName) %>%
  do(mean_and_category(.))

# Make try_byspmean into a species x trait matrix.

try_wide <- dcast(try_byspmean[,c(1,2,4)], AccSpeciesName ~ TraitName, fun.aggregate = mean)
try_wide <- try_wide[, apply(try_wide, 2, function(x) sum(!is.na(x))) > 10]

write.csv(try_wide, file.path(fp, 'try_fia.csv'), row.names = FALSE)

#################################################

# Try to run FD using the distance matrix generated

library(FD)

try_wide <- read.csv(file.path(fp, 'try_fia.csv'), stringsAsFactors = FALSE)

trydist <- gowdis(try_wide[, -1]) # Only two entries are NA which is pretty good. Can we resolve this?
data.frame(species = try_wide$AccSpeciesName, ntraits = apply(try_wide[,-1], 1, function(x) sum(!is.na(x))))

trydist2 <- gowdis(try_wide[!grepl('Cupressus', try_wide$AccSpeciesName), -1])


# Use only meaningful traits?
try_use <- try_wide[,c(3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 20, 21, 22, 24, 25)]
trydist <- gowdis(try_use)

ntraits <- data.frame(species = try_wide$AccSpeciesName, ntraits = apply(try_use, 1, function(x) sum(!is.na(x))))
ntraits %>% arrange(ntraits)


problemspp <- ntraits$species[ntraits$ntraits < 3] # We have to get rid of 15 species that don't have any traits. They seem rare though.

trydist2 <- gowdis(try_use[!try_wide$AccSpeciesName %in% problemspp, ])

# A good solution would be to just ignore the rare spp that don't have traits, and then not calculate FD for any of the plots that have at least 25% of those species abundance.

#################

# Load the fia matrix and try to figure out which species need traits.
fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'

load(file.path(fp, 'fiamatrices.RData'))

library(FD)

try_wide <- read.csv(file.path(fp, 'try_fia.csv'), stringsAsFactors = FALSE)


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

# Try to run FD on the fiaplotmat/tryspp combination.
try_use <- try_use[dimnames(fiaplotmat)[[2]], ]
trydist <- gowdis(try_use)
ntraits <- apply(try_use, 1, function(x) sum(!is.na(x)))

# Check which plots have greater than 25% species that we don't have traits for. Throw them out later.
problemspp <- names(ntraits)[ntraits < 3]
percent_problem <- apply(fiaplotmat, 1, function(x) sum(x[names(x) %in% problemspp])/sum(x))
# Only a few hundred plots have greater than 25% so we can just throw them out.

fiaplotmat_noproblem <- fiaplotmat[, !dimnames(fiaplotmat)[[2]] %in% problemspp]
try_noproblem <- try_use[!dimnames(try_use)[[1]] %in% problemspp, ]

# Code also breaks if a plot has zero trees in it.
zerorows <- apply(fiaplotmat_noproblem, 1, sum) == 0

fd_plot <- dbFD(x = try_noproblem, a = fiaplotmat_noproblem[!zerorows, ], w.abun = TRUE, corr = 'cailliez')

# Output of FD should have NAs for the zero-abundance communities, and the ones that were thrown out for having over 25% species with unknown traits.
n <- nrow(fiaplotmat)
fd_out <- data.frame(FRic = rep(NA, n), FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA)
fd_out$FRic[!zerorows] <- fd_plot$FRic
fd_out$FEve[!zerorows] <- fd_plot$FEve
fd_out$FDiv[!zerorows] <- fd_plot$FDiv
fd_out$FDis[!zerorows] <- fd_plot$FDis
fd_out$RaoQ[!zerorows] <- fd_plot$RaoQ

fd_out[percent_problem >= 0.25, ] <- NA

write.csv(fd_out, file.path(fp, 'fia_fd.csv'), row.names = FALSE)

#######
# Run FD at subplot level.

fiasubpmat_noproblem <- fiasubpmat[, !dimnames(fiasubpmat)[[2]] %in% problemspp]
zerorows_subplot <- apply(fiasubpmat_noproblem, 1, sum) == 0
percent_problem_subplot <- apply(fiasubpmat, 1, function(x) sum(x[names(x) %in% problemspp])/sum(x))

fd_subplot <- dbFD(x = try_noproblem, a = fiasubpmat_noproblem[!zerorows_subplot, ], w.abun = TRUE, corr = 'cailliez')

n <- nrow(fiasubpmat)
fd_outsubp <- data.frame(FRic = rep(NA, n), FEve = NA, FDiv = NA, FDis = NA, RaoQ = NA)
fd_outsubp$FRic[!zerorows_subplot] <- fd_subplot$FRic
fd_outsubp$FEve[!zerorows_subplot] <- fd_subplot$FEve
fd_outsubp$FDiv[!zerorows_subplot] <- fd_subplot$FDiv
fd_outsubp$FDis[!zerorows_subplot] <- fd_subplot$FDis
fd_outsubp$RaoQ[!zerorows_subplot] <- fd_subplot$RaoQ

fd_outsubp[percent_problem_subplot >= 0.25, ] <- NA
write.csv(fd_outsubp, file.path(fp, 'fia_fd_subplot.csv'), row.names = FALSE)
