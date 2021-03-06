# Compile spatial locations of trees in the try dataset with their trait information..
# These will be used for imputation.
# QDR 26 June 2017

####################################################################################

# TRY: generate data frame with spatial locations of all individuals from which traits were measured, with all their trait values and indications of whether traits are missing

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/'

try_all <- read.delim(file.path(fp, 'nasa_try.txt'),stringsAsFactors = FALSE, quote = '')

nmeas <- table(try_all$AccSpeciesName, try_all$TraitName)
unittable <- table(try_all$TraitName, try_all$OrigUnitStr)

# Create new column called flag and fill with NA's
try_all["flag"] <- NA

flagLatLong <- function(x)
{
  # If Longitude/Latitude are zeroes, then flag the row
  x$flag[((x$DataName == "Latitude" | x$DataName == "Longitude") &
            (x$OrigValueStr == '0'))] <- 'flag coordinates'
  return(x)
}

# Calls flagLatLong function on the data
try_all <- flagLatLong(try_all)

# Remove some of the columns that have irrelevant information.
# Numeric traits use StdValue instead of OrigValueStr
useStdValueIfNumeric <- function (x)
{
  x$OrigValueStr[!is.na.data.frame(as.numeric(x$OrigValueStr))] <- NA
  return(x)
}

try_all <- useStdValueIfNumeric(try_all)

#try_all <- try_all[,c('DatasetID','AccSpeciesName','ObservationID','TraitName','DataName','OrigValueStr','UnitName','StdValue', 'OrigUncertaintyStr','UncertaintyName')]

try_all$correct_value <- try_all$OrigValueStr # Writes character values to correct value column
try_all$correct_value[is.na(try_all$correct_value)] <- try_all$StdValue[is.na(try_all$correct_value)] # Writes Std Value to correct value

#######################################################

# Create one data frame with trait data rows only, and one with spatial data rows only. 
# Ignore the rest.

try_nometa <- subset(try_all, TraitName != '')
try_locations <- subset(try_all, DataName %in% c('Longitude','Latitude') & is.na(flag)) # removes flagged coordinates

# Reshape trait locations to only have the ID information and a separate column for latitude and longitude.

library(dplyr)

get_coords <- function(x) {
  lats <- x$StdValue[x$DataName == 'Latitude']
  lons <- x$StdValue[x$DataName == 'Longitude']
  lat_mode <- lats[which.max(table(lats))]
  lon_mode <- lons[which.max(table(lons))]
  data.frame(lat = lat_mode, lon = lon_mode)
}

try_locations_wide <- try_locations %>%
  group_by(AccSpeciesName, DatasetID, ObservationID) %>%
  do(get_coords(.))

# Get corrected individual measurements for each ObservationID.


# Figure out whether individual traits have more than one unit of measurement.
measByUnitTable <- table(try_nometa$TraitName, try_nometa$UnitName)
measByUnitTable[apply(measByUnitTable>0, 1, sum) > 1, ]

# Plant longevity has some blank units and some in years
longevitynounit <- subset(try_nometa, grepl('Plant lifespan',TraitName) & UnitName=='')
# Replace the names with two different values
try_nometa$TraitName[grepl('Plant lifespan',try_nometa$TraitName) & try_nometa$UnitName==''] <- 'Plant lifespan categorical'
try_nometa$TraitName[grepl('Plant lifespan',try_nometa$TraitName) & try_nometa$UnitName=='year'] <- 'Plant lifespan years'

# Seedbank longevity has some entries given as a percentage, some dimensionless, and some in years
seedbanknounit <- subset(try_nometa, grepl('\\(seedbank\\) longevity',TraitName) & UnitName=='dimensionless') # Categorical (e.g. transient, persistent)

try_nometa$TraitName[grepl('\\(seedbank\\) longevity',try_nometa$TraitName) & try_nometa$UnitName=='dimensionless'] <- 'Seedbank longevity categorical'

# Species by trait table for all the traits 
spByTraitTable <- table(try_nometa$AccSpeciesName, try_nometa$TraitName)
apply(spByTraitTable>0, 2, sum) # How many species do we have for each trait?

library(reshape2) # Package for reshaping data frames and matrices into long or wide format.

# Function to be run on each column to either return the mean or the first character result.
# Use most common character value rather than the first one.
mean_with_char <- function(x) {
  xnum <- as.numeric(x)
  xtable <- table(x, useNA = 'always')
  if (any(!is.na(xnum))) x <- as.character(mean(xnum, na.rm=TRUE)) else x <- names(xtable)[which.max(xtable)]
  return(x)
}

# Function to change columns in a df that only consist of numeric strings to numerics.
numstring2num <- function(x) {
  xnum <- as.numeric(x)
  if (any(!is.na(xnum))) xnum else x
}

# Cast (convert from long form to wide so that each trait has its own column and can be edited)

data_to_cast <- try_nometa[,c('ObservationID','DatasetID','AccSpeciesName','DataName','TraitName','correct_value')]

try_byobsmean <- dcast(data_to_cast, 
                       ObservationID+DatasetID+AccSpeciesName ~ TraitName, 
                       value.var='correct_value', 
                       fun.aggregate = mean_with_char)

# Convert all character columns in try_byobsmean that should be numeric to numeric.

try_byobsmean <- as.data.frame(lapply(try_byobsmean, numstring2num))

# Join trait and spatial info.

try_location_trait_byobs <- left_join(try_locations_wide, try_byobsmean)

write.csv(try_location_trait_byobs, file = 'X:/data/fia/try_location_trait_byobs.csv', row.names = FALSE)


######################################
# 30 June 2017: Generate different trait datasets with and without missing values to test imputation methods.

traits_to_use <- c(11, 42, 45, 58, 60)
traits_reduced <- try_location_trait_byobs[,c(1:5, traits_to_use)]
names(traits_reduced)[6:10] <- c('specific_leaf_area', 'plant_height', 'plant_lifespan', 'rooting_depth', 'seed_dry_mass')

# Use only individuals that have at least 1 of those traits measured.
userows <- apply(traits_reduced, 1, function(x) any(!is.na(x[6:10])))

traits_individual <- traits_reduced[userows,]

# Species means
traits_speciesmean <- traits_individual %>%
  group_by(AccSpeciesName) %>%
  summarize_at(4:10, mean, na.rm = TRUE)

write.csv(traits_individual, file = 'X:/jay/tree_data/traits_individual.csv', row.names=F)
write.csv(traits_speciesmean, file = 'X:/jay/tree_data/traits_speciesmean.csv', row.names=F)

# No missing values
traits_nomiss <- traits_speciesmean[,c(1:5,8)] %>% filter(complete.cases(.))
