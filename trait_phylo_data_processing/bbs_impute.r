# BBS numbers for generation time and for migrant status
# QDR/NASAbioXgeo/20 Jun 2018

# First, impute any missing longevity values, using the full trait dataset.
# Use phylogenetic imputation.

library(ape)
library(Rphylopars)

# Load trait data and separate into species IDs and trait data frames
bird_traits <- read.csv('/mnt/research/aquaxterra/DATA/raw_data/BBS/bird_traits/birdtraitmerged.csv', stringsAsFactors = FALSE)
bird_traits[bird_traits == -999] <- NA

use_traits <- c("Diet.Inv", "Diet.Vend", "Diet.Vect", "Diet.Vfish", "Diet.Vunk",
                "Diet.Scav", "Diet.Fruit", "Diet.Nect", "Diet.Seed", "Diet.PlantO",
                "ForStrat.watbelowsurf", "ForStrat.wataroundsurf", "ForStrat.ground",
                "ForStrat.understory", "ForStrat.midhigh", "ForStrat.canopy",
                "ForStrat.aerial", "PelagicSpecialist", "female_maturity_d",
                "litter_or_clutch_size_n", "litters_or_clutches_per_y", "adult_body_mass_g",
                "maximum_longevity_y", "birth_or_hatching_weight_g", "egg_mass_g",
                "incubation_d", "fledging_age_d", "longevity_y", "male_maturity_d"
)

bird_trait_ids <- bird_traits[,1:14]
bird_traits_use <- bird_traits[,use_traits]

# Load consensus tree for imputing
eric_cons_tree <- read.tree('/mnt/research/nasabio/data/bbs/ericson_cons.tre')

bird_trait_ids$Latin_Name_clean <- gsub(' ', '_', bird_trait_ids$Latin_Name_clean)

table(bird_trait_ids$Latin_Name_clean %in% eric_cons_tree$tip.label)

