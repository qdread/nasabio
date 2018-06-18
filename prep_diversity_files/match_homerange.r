# Match home range with species list

spp <- read.csv('~/GitHub/aquaxterra/data/specieslist.csv', stringsAsFactors = FALSE)
ranges <- read.csv('~/R/Tamburelloetal_HomeRangeDatabase.csv', stringsAsFactors = FALSE)

library(dplyr)

# Keep only birds, capitalize genus name, make sure species name is lowercase, and paste together genus and species to one string
# Also get rid of any extra white space
ranges <- ranges %>%
  filter(taxon == 'birds') %>%     
  mutate(genus = trimws(Hmisc::capitalize(genus)),
         species = trimws(tolower(species)),
         Latin_name = paste(genus, species))

# See which of the 140 bird species are not in our species list, using all possible synonyms
# 80 species don't match
no_match <- ranges$Latin_name[!ranges$Latin_name %in% spp$Latin_Name_clean &
                              !ranges$Latin_name %in% spp$Latin_Name_synonym &
                              !ranges$Latin_name %in% spp$Latin_Name_synonym2]

# Correct genus name typos
ranges$genus[grep('Laniuis', ranges$genus)] <- 'Lanius'
ranges$genus[grep('Typmanuchus', ranges$genus)] <- 'Tympanuchus'
ranges$genus[grep('Geothylpis|Geothlypis', ranges$genus)] <- 'Geothlypis'

# Recreate the Latin name column

ranges <- mutate(ranges, Latin_name = paste(genus, species))

no_match <- ranges$Latin_name[!ranges$Latin_name %in% spp$Latin_Name_clean &
                                !ranges$Latin_name %in% spp$Latin_Name_synonym &
                                !ranges$Latin_name %in% spp$Latin_Name_synonym2]
# Now down to 77 species that don't match

# The only ones I can see in here that are North American are Vireo, Setophaga, and Seiurus.
# I think their names are just spelled wrong in Tamburello database.

# Correct species name typos

ranges$Latin_name[grep('Vireo belli', ranges$Latin_name)] <- 'Vireo bellii'
ranges$Latin_name[grep('Setophaga kirtlandi', ranges$Latin_name)] <- 'Setophaga kirtlandii'
ranges$Latin_name[grep('Seiurus aurocapillus', ranges$Latin_name)] <- 'Seiurus aurocapilla'
ranges$Latin_name[grep('Vireo atricapillus', ranges$Latin_name)] <- 'Vireo atricapilla'

no_match <- ranges$Latin_name[!ranges$Latin_name %in% spp$Latin_Name_clean &
                                !ranges$Latin_name %in% spp$Latin_Name_synonym &
                                !ranges$Latin_name %in% spp$Latin_Name_synonym2]
# Now down to 73 species that don't match


# Create a column in the species list that has only the species names that match Tamburello
spp <- spp %>%
  mutate(Latin_name_for_joining_ranges = case_when(Latin_Name_clean %in% ranges$Latin_name ~ Latin_Name_clean,
                                                   Latin_Name_synonym %in% ranges$Latin_name ~ Latin_Name_synonym,
                                                   Latin_Name_synonym2 %in% ranges$Latin_name ~ Latin_Name_synonym2,
                                                   TRUE ~ Latin_Name_clean))

# Join species list with home range data by that column. For now just take a few columns from home range data

spp <- spp %>%
  left_join(ranges[,c('Latin_name', 'mean.hra.m2')], by = c(Latin_name_for_joining_ranges = 'Latin_name'))

table(!is.na(spp$mean.hra.m2)) # 67 matches
