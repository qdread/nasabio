# Test spatial imputation method using the cleaned up TRY data.
# QDR 5 July 2017

# Load TRY data that has valid coordinates, and remove some of the trait values at random.
try_geo <- read.csv('C:/Users/Q/google_drive/NASABiodiversityWG/Trait_Data/fia_try_04jul/try_trait_byobs_georeferenced.csv', stringsAsFactors = FALSE)

# Reduce to a few traits
library(dplyr)
try_geo_sub <- try_geo %>%
  select(AccSpeciesName, ObservationID, lat, lon, Leaf.area.per.leaf.dry.mass..specific.leaf.area..SLA., Plant.height, Seed.dry.mass)

names(try_geo_sub) <- c('species','individual','lat','lon','SLA','height','seedmass')
try_geo_complete <- filter(try_geo_sub, complete.cases(try_geo_sub))
try_geo_atleast2 <- filter(try_geo_sub, (as.numeric(!is.na(SLA)) + as.numeric(!is.na(height)) + as.numeric(!is.na(seedmass))) >= 2)
# This might not be the best dataset but it will have to do.

# Load spatial imputation
library(yaImpute)
