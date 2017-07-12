# environmental predictors, phylogenetic VCV matrix, trait data, and spatial locations
library(dplyr)

# For now, create environmental predictors and traits from species means (mean locations)

# Use Stevens traits, and get median locations of each species from FIA dataset.
# Use the environmental predictors from FIA.

trait_training <- read.csv('X:/jay/tree_data/stevenstraits/trait_stevens_training.csv', stringsAsFactors = FALSE)

# Get median locations for each species from TRY
# All covariates for those locations have already been extracted from different GIS layers.
try_covariates <- read.csv('X:/data/fia/imputation_covariates.csv', stringsAsFactors = FALSE)

# Get only a subset of the covariates: elevation, minimum temperature, temperature seasonality, minimum precipitation, precipitation seasonality, lai, npp
try_covariates_reduced <- try_covariates %>%
  select(AccSpeciesName, lat_aea, lon_aea, elevation, bio6, bio4, bio14, bio15, lai, npp) %>%
  rename(temp_min_monthly = bio6, temp_seasonality = bio4, precip_min_monthly = bio14, precip_seasonality = bio15) %>%
  filter(complete.cases(.))


library(ape)
allfia_phylogeny <- read.tree('X:/data/fia/tree_all_final_031716.txt')

# Pare down the trait list and phylogeny to species that we have valid locations for. (64 spp)

species_list <- Reduce('intersect', list(try_covariates_reduced$AccSpeciesName, allfia_phylogeny$tip.label, trait_training$Scientific_Name))

species_covariates <- subset(try_covariates_reduced, AccSpeciesName %in% species_list)
species_phylogeny <- drop.tip(allfia_phylogeny, tip = allfia_phylogeny$tip.label[!allfia_phylogeny$tip.label %in% species_list])
species_traits <- subset(trait_training, Scientific_Name %in% species_list)

# Sort them all so that they are in the same order
species_covariates <- species_covariates %>%
  rename(species = AccSpeciesName) %>%
  arrange(species)

species_traits <- species_traits %>%
  rename(species = Scientific_Name) %>%
  arrange(species)