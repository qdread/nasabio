# Load plant phylogeny and extract species in FIA to calculate phylogenetic diversity of FIA plots.

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'
fiapnw <- read.csv(file.path(fp, 'finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)

library(ape)

planttree <- read.tree(file.path(fp, 'embryophyta.lter.elbowsremoved.tre')) # Has all spp in LTER so should have all species in FIA plots, except maybe a few oddballs.

# Get FIA species code to scientific name lookup table, and find which ones are in the tree.
fiataxa <- read.csv('specieslists/fia_taxon_lookuptable.csv', stringsAsFactors = FALSE)
pnw_codes <- unique(fiapnw$SPCD)
all(pnw_codes %in% fiataxa$FIA.Code) #yes

# Scientific names of FIA species in PNW plots.
pnw_species <- fiataxa[match(pnw_codes, fiataxa$FIA.Code), c('Genus','Species')]
pnw_scinames <- paste(pnw_species$Genus, pnw_species$Species, sep = '_')
pnw_scinames <- gsub(' ', '', pnw_scinames) #Remove extraneous spaces

phymatch <- pnw_scinames %in% planttree$tip.label
pnw_scinames[!phymatch]

# Check whether at least all the genera are represented.
pnw_genera <- unique(gsub(' ', '', pnw_species$Genus))
tree_genera <- sapply(strsplit(planttree$tip.label, '_'), '[', 1)
pnw_genera[!pnw_genera %in% tree_genera]
# Actually missing quite a few from this phylogeny. Might be better to get another phylogeny.

# Got phylogeny from Smith et al. 2011, Am J Bot: http://datadryad.org/resource/doi:10.5061/dryad.8790

planttree2 <- read.tree(file.path(fp, 'final_tree.tre'))
phymatch2 <- pnw_scinames %in% planttree2$tip.label
pnw_scinames[!phymatch2]

tree_genera2 <- sapply(strsplit(planttree2$tip.label, '_'), '[', 1)
pnw_genera[!pnw_genera %in% tree_genera2]
# All the genera are represented in this new phylogeny so that's good.

# Add non-matching names to the tree in the correct genus.
# This approximation should not be a large problem.

######
# example for mammals:

library(phytools)
for (i in spp_not_in_tree) mammaltol <- add.species.to.genus(tree = mammaltol, species = i, where = 'root')
#for (i in spp_not_in_tree[10:13]) mammaltol <- add.species.to.genus(tree = mammaltol, species = i, where = 'random')

