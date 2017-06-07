# Output TRY species IDs for all FIA species for new TRY request
# 7 June 2017

# Load FIA species list (binomials)
fia_spp <- read.csv('specieslists/fia_taxon_lookuptable.csv', stringsAsFactors = FALSE)

# Select only "core" species (exclude the Puerto Rican species)
fia_spp <- subset(fia_spp, Core != '')

# Load TRY species list
try_spp <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/tryspp.csv', stringsAsFactors = FALSE)

# Match FIA species list with TRY species list
# Need to get rid of spaces in the FIA species names
# Also fix two typos in TRY
try_spp$AccSpeciesName[try_spp$AccSpeciesName == 'Quercus margaretta'] <- 'Quercus margarettae'
try_spp$AccSpeciesName[try_spp$AccSpeciesName == 'Gymnocladus dioica'] <- 'Gymnocladus dioicus'
try_spp$AccSpeciesName[try_spp$AccSpeciesName == 'Casuarina pauper'] <- 'Casuarina lepidophloia'


library(dplyr)
fia_spp <- fia_spp %>%
  mutate(Genus = gsub(" ", "", Genus, fixed = TRUE),
         Species = gsub(" ", "", Species, fixed = TRUE),
         AccSpeciesName = paste(Genus, Species)) %>%
  left_join(try_spp)

summary(fia_spp) # All but 10 rows have a valid ID
fia_spp$AccSpeciesName[is.na(fia_spp$AccSpeciesID)] # They are all either subspecies or unknown.

# Display species IDs
paste(na.omit(fia_spp$AccSpeciesID), collapse = ',')

# Save to CSV
write.csv(fia_spp, 'specieslists/fia_try_spp.csv', row.names = FALSE)
