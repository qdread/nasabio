try_all <- read.delim('C:/Users/Q/Dropbox/projects/nasabiodiv/fia_try_04jul2017.txt',stringsAsFactors = FALSE, quote = '')

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
  lat_mode <- lats[which.max(table(lats, useNA = 'always'))]
  lon_mode <- lons[which.max(table(lons, useNA = 'always'))]
  data.frame(lat = lat_mode, lon = lon_mode)
}

try_locations_wide <- try_locations %>%
  group_by(AccSpeciesName, DatasetID, ObservationID) %>%
  do(get_coords(.))

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

# LDMC has some with no unit
ldmcnounit <- subset(try_nometa, grepl('LDMC',TraitName) & UnitName=='')
# They are just missing values
try_nometa <- subset(try_nometa, !(grepl('LDMC',TraitName) & UnitName==''))

# Sclerophylly has some with no unit. Only a few are in N/mm
scleronounit <- subset(try_nometa, grepl('sclerophylly',TraitName) & UnitName=='')
# They are a variety of things.
try_nometa$TraitName[grepl('sclerophylly',try_nometa$TraitName) & try_nometa$UnitName==''] <- 'Sclerophylly categorical'

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

# Get a small subset that is actually good data.
try_testdata <- try_location_trait_byobs %>%
  rename(Specific.leaf.area = Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA.,
         Leaf.Cmass = Leaf.carbon..C..content.per.leaf.dry.mass,
         Leaf.dry.matter.content = Leaf.dry.mass.per.leaf.fresh.mass..Leaf.dry.matter.content..LDMC.,
         Leaf.Nmass = Leaf.nitrogen..N..content.per.leaf.dry.mass,
         Leaf.Pmass = Leaf.phosphorus..P..content.per.leaf.dry.mass,
         Specific.stem.density = Stem.dry.mass.per.stem.fresh.volume..stem.specific.density..SSD..wood.density.,
         Stomatal.conductance = Stomata.conductance.per.leaf.area) %>%
  select(AccSpeciesName, DatasetID, ObservationID, lat, lon, Specific.leaf.area, Leaf.thickness, Leaf.Cmass, Leaf.Nmass, Leaf.Pmass, Leaf.dry.matter.content, Plant.height, Seed.dry.mass, Specific.stem.density, Stomatal.conductance)

try_smalldataset <- try_testdata %>%
  ungroup %>%
  select(AccSpeciesName, lat, lon, Specific.leaf.area, Leaf.Nmass, Plant.height) %>%
  filter(complete.cases(.))

# Centroid of coords
get_centroid <- function(x) {
  coords <- cbind(x$lat, x$lon)
  uniquecoords <- unique(coords)
  data.frame(lat = mean(uniquecoords[,1], na.rm=T), lon = mean(uniquecoords[,2], na.rm=T))
}

# Get median locations for each species (leave out non-US locations)
try_spp_loc <- try_locations_wide %>%
  ungroup %>%
  filter(lon < -50, lat > 0) %>%
  group_by(AccSpeciesName) %>%
  do(get_centroid(.)) %>%
  ungroup %>%
  mutate(AccSpeciesName = gsub('\\ ', '_', AccSpeciesName))

write.csv(try_spp_loc, file = 'X:/data/fia/tree_locations.csv', row.names = FALSE)


# Get the unit names for each trait ---------------------------------------

try_unitnames <- try_nometa %>%
  group_by(TraitName) %>%
  summarize(unit = unique(UnitName))

write.csv(try_unitnames, file = 'C:/Users/Q/google_drive/NASABiodiversityWG/Trait_Data/try_units.csv', row.names = FALSE)
