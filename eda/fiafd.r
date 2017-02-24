# Functional diversity of FIA
fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'

# First get list of species IDs for TRY.
tryspp <- read.csv(file.path(fp, 'tryspp.csv'), stringsAsFactors = FALSE)

trymatch <- pnw_species$sciname %in% tryspp$AccSpeciesName
pnw_species$sciname[!trymatch] # All match except for the "unknowns"

tryids <- tryspp$AccSpeciesID[match(pnw_species$sciname, tryspp$AccSpeciesName)]
write.table(t(na.omit(tryids)), sep=',')

# Load the try species list.
try_all <- read.delim(file.path(fp, 'nasa_try.txt'), stringsAsFactors = FALSE)

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
seedbankpct <- subset(try_nometa, grepl('\\(seedbank\\) longevity',TraitName) & UnitName=='%') # Fraction of plots with persistent seeds

try_nometa$TraitName[grepl('\\(seedbank\\) longevity',try_nometa$TraitName) & try_nometa$UnitName=='dimensionless'] <- 'Seedbank longevity categorical'
try_nometa$TraitName[grepl('\\(seedbank\\) longevity',try_nometa$TraitName) & try_nometa$UnitName=='%'] <- 'Seedbank persistence percent'

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


library(readr)

try_byobsmean <- type_convert(try_byobsmean)
