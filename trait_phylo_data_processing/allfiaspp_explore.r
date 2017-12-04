# Script to explore FIA data for the whole entire USA.
# Run on cluster because it's a big file (900megs)

fiaall <- read.csv('finley_trees_continental_US_most_recent_evaluations_nov8_2017.csv', stringsAsFactors = FALSE)
fiaplotids <- read.csv('finley_unique_plt_cn_from_tree_dataset_nov8_2017.csv', stringsAsFactors = FALSE)
# Note: just loading fiaall into ram takes up ~4gb ram. Probably a good idea to split it up at some point. ~4 million rows, 29 columns.

fia_spp <- table(fiaall$SPCD)
length(fia_spp) # 383 species.

length(unique(fiaall$PLT_CN)) # 135174 plots.

# Load the coordinates
fiaallcoords <- read.csv('~/data/allfia.csv')
dim(fiaallcoords)
length(unique(fiaallcoords$CN)) # Both of these agree, 135174 plots. PNW is 1/6 of these, so I need to calc things for 5x as many plots as I've done before.

# Check the species codes against the lookup table.
# Save the list of species codes to do this locally.
allfia_spcodes <- unique(fiaall$SPCD)
write.csv(data.frame(SPCD = allfia_spcodes), file = 'spcds.csv', row.names = FALSE)

###################
# starting here, run locally from github working directory:

allfia_spcodes <- read.csv('specieslists/spcds.csv')
fiataxa <- read.csv('specieslists/fia_taxon_lookuptable.csv', stringsAsFactors = FALSE)
all(allfia_spcodes$SPCD %in% fiataxa$FIA.Code) #yes

# Browse the species list.
with(fiataxa, paste(Genus, Species)[match(allfia_spcodes$SPCD, FIA.Code)])
# This contains problematic entries: var and ssp, Genus with only sp given, Family Arecaceae, Tree evergreen, Tree broadleaf, Tree unknown

# See what we need to do to make sure all species are included in the trait db and the phylogeny.

library(ape)
library(dplyr)
fullphylo <- read.tree('C:/Users/Q/Dropbox/projects/nasabiodiv/allfiaphylogeny/tree_all_final_031716.txt')

fiataxa <- fiataxa %>%
  mutate(Genus = gsub('\\ ', '', Genus),
         Species = gsub('\\ ', '', Species),
         binomial = paste(Genus, Species, sep = '_'))

fiataxa_inplots <- fiataxa$binomial[match(allfia_spcodes$SPCD, fiataxa$FIA.Code)]
phylotaxa <- fullphylo$tip.label

table(fiataxa_inplots %in% phylotaxa)
table(phylotaxa %in% fiataxa_inplots)
# 337 overlap, but 46 of the taxa in the FIA plots aren't in the phylogeny.

fiataxa_inplots[!(fiataxa_inplots %in% phylotaxa)]
phylotaxa[!(phylotaxa %in% fiataxa_inplots)]

grep('^([^_]*_){2}[^_]*$', phylotaxa, value = TRUE) # Find all with exactly 2 underscores, see https://stackoverflow.com/questions/863125/regular-expression-to-count-number-of-commas-in-a-string.

# Issues that can easily be fixed (9):
# Carya carolinae-septentrionalis, Crataegus crus-galli: hyphen in fia, underscore in phylogeny
# Quercus margarettae: misspelled in phylogeny
# Tilia americana caroliniana and heterophylla: Replace the var. with an underscore.
# Castanea pumila ozarkensis: Replace the var. with an underscore.
# Populus balsamifera trichocarpa: Replace the var. with an underscore, correct capitalization
# Populus deltoides monilifera: Replace the var. with an underscore, correct capitalization
# Sapindus saponaria drummondii: Replace the var. with an underscore

# Trees classified at a level below species that we will need to ignore the below-species level (8):

# Juniperus virginiana silicicola, Quercus sinuata sinuata, Sideroxylon lanuginosum lanuginosum, Abies lasiocarpa arizonica, Chrysolepis chrysophylla chrysophylla, Pinus monophylla fallax, Aesculus glabra arguta, Celtis laevigata reticulata

# Trees not identified to species that we will either have to throw out or somehow use the coarser genus or family level ID (23)
# Amelanchier sp, Prunus sp, Quercus sp, Fraxinus sp, Diospyros sp, Family Arecaceae, Ulmus sp, Carya sp, Morus sp, Larix sp, Tamarix sp, Populus sp, Magnolia sp, Celtis sp, Sorbus sp, Pinus sp, Juniperus sp, Acer sp, Picea sp, Abies sp, Aesculus sp, Prosopis sp, Betula sp

### example of how to deal with the ones above. This probably could be modified to add Arecaceae sp. too.
# library(phytools)
# fullphylo <- add.species.to.genus(tree = fullphylo, species = 'Fraxinus_sp', where = 'root')

# Unknown trees that will have to be thrown out (3): Tree unknown, Tree broadleaf, Tree evergreen

# Trees that just aren't in the phylogeny (non-native species) (2): Prunus cerasus (might be able to assign to P. avium), Ginkgo biloba (nah fam)

# Trees whose species was lumped with another species and can be reassigned to that species (1): Abies shastensis (now under A. magnifica)

############################################

# Correct taxa names to match phylogeny.
# 1. correct capitalization and replace the subspecies and variety with underscores, replace the hyphens with underscores too
fiataxa <- fiataxa %>%
  mutate(Species = tolower(Species),
         Species = gsub('var\\.', '_', Species),
         Species = gsub('ssp\\.', '_', Species),
         Species = gsub('\\-', '_', Species),
         binomial = paste(Genus, Species, sep = '_'))

fiataxa_inplots <- fiataxa$binomial[match(allfia_spcodes$SPCD, fiataxa$FIA.Code)]

table(fiataxa_inplots %in% phylotaxa)
table(phylotaxa %in% fiataxa_inplots)

# 2. Get rid of the below-species level classifications that are not in the phylogeny.
spp_to_fix <- fiataxa_inplots[!fiataxa_inplots %in% phylotaxa]

badsubspecies <- grep('^([^_]*_){2}[^_]*$', spp_to_fix, value = TRUE)

# remove the 2nd underscore and everything following it.
fixedsubspecies <- sapply(strsplit(badsubspecies, '_'), function(x) paste(x[1:2], collapse = '_'))

fiataxa <- mutate(fiataxa, binomial_forphylo = binomial)
fiataxa$binomial_forphylo[match(badsubspecies, fiataxa$binomial_forphylo)] <- fixedsubspecies

fiataxa_inplots <- fiataxa$binomial_forphylo[match(allfia_spcodes$SPCD, fiataxa$FIA.Code)]

table(fiataxa_inplots %in% phylotaxa)

# 3. Add genera
genera_to_add <- grep('_sp$', spp_to_fix, value = TRUE)
library(phytools)
for (i in genera_to_add) {
  fullphylo <- add.species.to.genus(tree = fullphylo, species = i, where = 'root')
}

# Check whether this was successful
genera_to_add %in% fullphylo$tip.label # One didn't work
genera_to_add[!genera_to_add %in% fullphylo$tip.label] # Tamarix just isn't in the phylogeny either. Probably throw it out. It's only represented by invasives and probably is almost never in "forested" plots.

phylotaxa <- fullphylo$tip.label
table(fiataxa_inplots %in% phylotaxa)

spp_to_fix <- fiataxa_inplots[!fiataxa_inplots %in% phylotaxa] 

# Correct typo in phylogeny
fullphylo$tip.label[fullphylo$tip.label == 'Quercus_margarettiae'] <- 'Quercus_margarettae'

# Assign Prunus cerasus to Prunus avium, and assign Abies shastensis to Abies magnifica.
fiataxa$binomial_forphylo[fiataxa$binomial_forphylo %in% c('Abies_shastensis', 'Prunus_cerasus')] <- c('Abies_magnifica', 'Prunus_avium')

phylotaxa <- fullphylo$tip.label
fiataxa_inplots <- fiataxa$binomial_forphylo[match(allfia_spcodes$SPCD, fiataxa$FIA.Code)]

table(fiataxa_inplots %in% phylotaxa)

spp_not_in_phylo <- fiataxa_inplots[!fiataxa_inplots %in% phylotaxa] 

############################################

# Traits.
# The imputation was only done for the Pacific Northwest trees. I have to redo it for all the trees.
