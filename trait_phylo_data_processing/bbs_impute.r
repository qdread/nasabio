# BBS phylogenetic imputation of all traits
# QDR/NASAbioXgeo/20 Jun 2018

# Use phylogenetic imputation.

# Last modified on 28 Jun 2018: Finished replacing all traits of subspecies and hybrids with appropriate values.

library(ape)
library(Rphylopars)
library(dplyr)

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

tree_spp <- gsub('_', ' ', eric_cons_tree$tip.label) # Version of species in phylogeny with underscores removed

# We already have up to 3 possible synonyms for each bird species name, based on previous work I did. 
# Using all possible synonyms, let's see how many of them are in the phylogeny.

matchids <- with(bird_trait_ids,
     Latin_Name_clean %in% tree_spp | Latin_Name_synonym %in% tree_spp | Latin_Name_synonym2 %in% tree_spp)
table(matchids) # 94 species still don't match.

# Use the original latin names to use for matching the other 94.
nomatchnames <- bird_trait_ids$Latin_Name[!matchids] # Most are subspecies, hybrids, and unknowns. I already came up with a classification for them
table(bird_trait_ids$Type[!matchids])

# for phylogenetic purposes, do the following:
# Before imputing:
# 1. try to correct the ones that are "good" species
# 2. get rid of the ones that can only be identified to family level, and the species that are just plain not in it. There is no hope for them.
# 3. add genera to the phylogeny at the root so that "genus sp." can be identified.
# After imputing:
# 4. assign subspecies to the "parent" species, after doing the imputation
# 5. give hybrids the mean of the parent species traits, after doing the imputation
# 6. give the species that could be narrowed down to two species the mean of those two species, after doing the imputation

bird_trait_ids$Latin_Name[!matchids & bird_trait_ids$Type=='species']


#### create list of names for phylogeny
species_names_for_phylo <- case_when(
  bird_trait_ids$Latin_Name_synonym2 %in% tree_spp ~ bird_trait_ids$Latin_Name_synonym2,
  bird_trait_ids$Latin_Name_synonym %in% tree_spp ~ bird_trait_ids$Latin_Name_synonym,
  bird_trait_ids$Latin_Name_clean %in% tree_spp ~ bird_trait_ids$Latin_Name_clean,
  TRUE ~ as.character(NA)
)

# Step 1 from above list:
#### Correct species that have other synonyms that I didn't find before
new_ids <- rep(NA, nrow(bird_trait_ids))
new_ids[which(bird_trait_ids$Latin_Name_clean == 'Artemisiospiza nevadensis')] <- 'Amphispiza belli'
new_ids[which(bird_trait_ids$Latin_Name_clean == 'Troglodytes hiemalis')] <- 'Troglodytes troglodytes'
new_ids[which(bird_trait_ids$Latin_Name_clean == 'Thryophilus sinaloa')] <- 'Thryothorus sinaloa'
new_ids[which(bird_trait_ids$Latin_Name_clean == 'Antrostomus arizonae')] <- 'Caprimulgus vociferus'

species_names_for_phylo[!is.na(new_ids)] <- new_ids[!is.na(new_ids)]

# Step 2 from above list just means those ones won't be included in the imputation anyway.

# Step 3 from above list:
# Add "xxx sp." to phylogeny at the root of genus "xxx" using phytools functions
library(phytools)

genera_to_add <- bird_trait_ids$Latin_Name[bird_trait_ids$Type == 'unknown_genus']
# Some of these are only to the family level. Get only the ones that are to the genus level. (input manually)
genera_to_add <- c("Gavia sp.", "Phalacrocorax sp.", "Accipiter sp.", "Buteo sp.", "Calidris sp.", "Sphyrapicus sp.", "Empidonax sp.", "Poecile sp.")
genera_to_add_formatted <- gsub(' sp.', '_sp', genera_to_add, fixed = TRUE) # Use underscore formatting

# Make tree meet the assumptions needed to do the imputation
# Evolutionary biologists might not like these quick and dirty corrections, but they're fine for our purposes.
eric_cons_tree$root.edge <- 0
eric_cons_addedgenera <- force.ultrametric(eric_cons_tree, method = 'nnls')

for (i in genera_to_add_formatted) {
  eric_cons_addedgenera <- add.species.to.genus(tree = eric_cons_addedgenera, species = i, where = 'root')
}

# Check whether this was successful
genera_to_add_formatted %in% eric_cons_addedgenera$tip.label # One didn't work. This is because the old name Parus for Poecile is in the tree.

# Add Parus_sp
eric_cons_addedgenera <- add.species.to.genus(tree = eric_cons_addedgenera, species = 'Parus_sp', where = 'root')
species_names_for_phylo[which(bird_trait_ids$Latin_Name == 'Poecile sp.')] <- 'Parus sp.'

# Add the genus names to the phylogeny species name list
for (i in genera_to_add[1:7]) {
  species_names_for_phylo[which(bird_trait_ids$Latin_Name == i)] <- i
}

#########################################################
# Do the imputation
# #######################################################

# Format species names the same as phylogeny (get rid of periods and replace underscores with spaces)
species_names_for_phylo_formatted <- gsub('.', '', species_names_for_phylo, fixed = TRUE)
species_names_for_phylo_formatted <- gsub(' ', '_', species_names_for_phylo_formatted, fixed = TRUE)

library(Rphylopars)
# Get only the tree with the species in the trait dataset. (666 species, oooo)
imputation_tree <- drop.tip(eric_cons_addedgenera, tip = eric_cons_addedgenera$tip.label[!eric_cons_addedgenera$tip.label %in% species_names_for_phylo_formatted])
imputation_tree$root.edge <- 0

# Get only the traits for the species in the imputation tree
imputation_traits <- bird_traits_use[match(imputation_tree$tip.label, species_names_for_phylo_formatted), ]

# Here we actually do the imputation. It takes a few minutes.
# Log-transform all mass or weight variables before running imputation.
# This uses all columns for imputation, including those with no missing values (helps with covariance)
set.seed(5555)
phyimp <- phylopars(trait_data = data.frame(species = imputation_tree$tip.label, imputation_traits) %>%
                      mutate_at(vars(matches('mass|weight')), funs(log(.))),
                    tree = imputation_tree,
                    model = 'OU')
apply(phyimp$anc_recon, 2, min) # check for validity

# Back-transform and combine everything together.
traits_imputed <- phyimp$anc_recon[1:nrow(imputation_traits),] %>%
  as.data.frame %>%
  mutate_at(vars(matches('mass|weight')), funs(exp(.)))

dimnames(traits_imputed)[[1]] <- dimnames(phyimp$anc_recon)[[1]][1:nrow(traits_imputed)]

# Make a full data frame with both the imputed values and the ones we still need to replace
# Only replace the rows where at least 1 value was imputed
traits_imputed_subset <- traits_imputed[apply(imputation_traits, 1, function(x) any(is.na(x))), ]
traits_imputed_full <- bird_traits_use
for (i in 1:nrow(traits_imputed_full)) {
  if(species_names_for_phylo_formatted[i] %in% dimnames(traits_imputed_subset)[[1]]) {
    traits_imputed_full[i, ] <- traits_imputed_subset[dimnames(traits_imputed_subset)[[1]] == species_names_for_phylo_formatted[i], ]
  }
}

# Post processing (steps 4-6 of the above list)

library(purrr)
# Step 4: Give subspecies the values of their parent species 
# Make two corrections to this list for Setophaga and Rallus

subspecies <- bird_trait_ids %>%
  filter(Type == 'subspecies') %>%
  mutate(name = paste(Genus, map_chr(strsplit(Species, ' '), 1))) %>%
  mutate(name = case_when(grepl('Setophaga', name) ~ 'Dendroica coronata',
                          grepl('Rallus', name) ~ 'Rallus longirostris',
                          TRUE ~ name))

# Replace the subspecies traits
traits_imputed_full[which(bird_trait_ids$Type == 'subspecies'), ] <- traits_imputed_full[match(subspecies$name, species_names_for_phylo), ]


# Step 5: Give hybrids the mean values of their parent species
# Use AOU number ID to match
hybrids <- bird_trait_ids %>%
  filter(Type == 'hybrid')
# Retrieve AOUs
hybrid_ids <- map(hybrids$AOU_list, strsplit, split = ' ') %>%
  map(1) %>%
  map(as.numeric)

# Get the mean values for each of those multiple AOUs
hybrid_traits <- map_dfr(hybrid_ids, function(ids) traits_imputed_full[match(ids, bird_trait_ids$AOU), ] %>% summarize_all(mean))

# Replace the hybrid traits
traits_imputed_full[which(bird_trait_ids$Type == 'hybrid'), ] <- hybrid_traits

# Step 6: Give individuals that could only be narrowed down to 2 or 3 species the mean values of those 2 or 3 species
# Again use AOU number ID to match
ambiguous <- bird_trait_ids %>%
  filter(Type == 'unknown_sp')
# Retrieve AOUs
ambiguous_ids <- map(ambiguous$AOU_list, strsplit, split = ' ') %>%
  map(1) %>%
  map(as.numeric)

# Get the mean values for each of those multiple AOUs
ambiguous_traits <- map_dfr(ambiguous_ids, function(ids) traits_imputed_full[which(bird_trait_ids$AOU %in% ids), ] %>% summarize_all(mean, na.rm = TRUE))

# Replace the hybrid traits
traits_imputed_full[which(bird_trait_ids$Type == 'unknown_sp'), ] <- ambiguous_traits

# See where any missing values are.
bird_trait_ids[!complete.cases(traits_imputed_full), ]
# The only ones left are three that are only identified to family!!! We are good!

###########################
# Write imputed trait data to file.

traits_imputed_withids <- data.frame(AOU =  bird_trait_ids$AOU,
                                     Latin_Name_from_traitdb = bird_trait_ids$Latin_Name,
                                     Latin_Name_clean = bird_trait_ids$Latin_Name_clean,
                                     Latin_Name_phylogeny = species_names_for_phylo,
                                     Family = bird_trait_ids$Family,
                                     Genus = bird_trait_ids$Genus,
                                     Species = bird_trait_ids$Species,
                                     English_Common_Name = bird_trait_ids$English_Common_Name,
                                     Type = bird_trait_ids$Type,
                                     traits_imputed_full)
write.csv(traits_imputed_withids, file = '/mnt/research/nasabio/data/bbs/bbs_traits_imputed.csv', row.names = FALSE)