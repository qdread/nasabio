# Stevens trait data and Potter phylogeny
# FIA species
# 30 June 2017

library(XLConnect)
library(dplyr)
trait_stevens <- readWorksheetFromFile('C:/Users/Q/google_drive/NASABiodiversityWG/Trait_Data/Traits_Stevens_FIA.xlsx', sheet = 'Master')

# Convert to numerics where needed.
numeric_cols <- c(2,5:ncol(trait_stevens))
trait_stevens[, numeric_cols] <- lapply(trait_stevens[,numeric_cols], as.numeric)

# Get rid of subspecies and genus-level trait values
trait_stevens <- filter(trait_stevens, Scientific_Name != 'Tree_unknown', Subsp %in% c(0,0.5), Generic != 1)



# Phylogeny linked with these
library(ape)
allfia_phylogeny <- read.tree('C:/Users/Q/Dropbox/projects/nasabiodiv/allfiaphylogeny/tree_all_final_031716.txt')

# How many of Stevens scientific names are in there.
spmatch <- trait_stevens$Scientific_Name %in% allfia_phylogeny$tip.label
table(spmatch)
trait_stevens$Scientific_Name[!spmatch]

# Correct names that are typos or obsolete names in trait_stevens
name_correction <- c('Chamaecyparis_lawsonia' = 'Chamaecyparis_lawsoniana',
                     'Aesculus_octandra' = 'Aesculus_flava',
                     "Carya_illinoensis" = "Carya_illinoinensis" ,
                     'Carya_tomentosa' = 'Carya_alba',
                     'Populus_trichocarpa' = "Populus_balsamifera_trichocarpa",
                     "Quercus_chrysolepsis" = 'Quercus_chrysolepis',
                     'Quercus_nutallii' = 'Quercus_texana',
                     'Quercus_wislizenii' = 'Quercus_wislizeni',
                     'Tilia_heterophylla' = 'Tilia_americana_heterophylla',
                     'Ulmus_pumilia' = 'Ulmus_pumila')

trait_stevens$Scientific_Name[match(names(name_correction), trait_stevens$Scientific_Name)] <- name_correction

# 3 still don't match
spmatch <- trait_stevens$Scientific_Name %in% allfia_phylogeny$tip.label
table(spmatch)
trait_stevens <- filter(trait_stevens, spmatch)

trait_stevens <- trait_stevens[,-c(1,2,3,5:12)]
trait_stevens_continuous <- trait_stevens[,c(1,2,3,6:19)]
trait_stevens_categorical <- trait_stevens[,c(1,4,5,20,21)]

# Small training dataset with no missing values. We can place missing values in there at random.
training_traits <- c('Bark.thickness', 'Wood.density', 'SLA', 'Plant.height', 'Plant.lifespan', 'Seed.dry.mass')
trait_stevens_training <- trait_stevens[, c('Scientific_Name',training_traits)] %>%
  filter(complete.cases(.))

# Save them all
write.csv(trait_stevens_training, file = 'X:/jay/tree_data/stevenstraits/trait_stevens_training.csv', row.names = FALSE)
write.csv(trait_stevens_categorical, file = 'X:/jay/tree_data/stevenstraits/trait_stevens_categorical.csv', row.names = FALSE)
write.csv(trait_stevens_continuous, file = 'X:/jay/tree_data/stevenstraits/trait_stevens_continuous.csv', row.names = FALSE)
