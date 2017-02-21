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

# Correction for subspecies that are not in the tree.
pnw_species$sciname <- pnw_scinames
pnw_species$sciname[pnw_species$sciname == 'Chrysolepis_chrysophyllavar.chrysophylla'] <- 'Chrysolepis_chrysophylla'
pnw_species$sciname[pnw_species$sciname == 'Populus_balsamiferassp.Trichocarpa'] <- 'Populus_trichocarpa'

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

spp_not_in_tree <- sort(pnw_species$sciname[!phymatch2])

# Don't use the totally unidentified species
spp_not_in_tree <- spp_not_in_tree[!grepl('Tree', spp_not_in_tree)]

######
# Adding species to planttree2 does not work because it does not have any branch lengths.

library(phytools)
for (i in spp_not_in_tree) planttree2 <- add.species.to.genus(tree = planttree2, species = i, where = 'random')

######
# Format species list for phylocom/phylomatic.
# Add families, and get rid of the unidentified species/subspecies

pnw_species <- fiataxa[match(pnw_codes, fiataxa$FIA.Code), c('Genus','Species')]
pnw_scinames <- paste(pnw_species$Genus, pnw_species$Species, sep = '_')
pnw_scinames <- gsub(' ', '', pnw_scinames) #Remove extraneous spaces

# Correction for subspecies that are not in the tree.
pnw_species$sciname <- pnw_scinames
pnw_species$sciname[pnw_species$sciname == 'Chrysolepis_chrysophyllavar.chrysophylla'] <- 'Chrysolepis_chrysophylla'
pnw_species$sciname[pnw_species$sciname == 'Populus_balsamiferassp.Trichocarpa'] <- 'Populus_trichocarpa'

pnw_species <- subset(pnw_species, !grepl('Tree', sciname)) # unidentified species
pnw_species$sciname <- gsub('_', ' ', pnw_species$sciname)
pnw_species <- pnw_species[order(pnw_species$sciname), ]

family_table <- c(Abies='Pinaceae', Acer='Sapindaceae', Aesculus='Sapindaceae', Ailanthus='Simaroubaceae', Alnus = 'Betulaceae', Arbutus = 'Ericaceae', Betula = 'Betulaceae', Calocedrus = 'Cupressaceae', Cercocarpus = 'Rosaceae', Chamaecyparis = 'Cupressaceae', Chrysolepis = 'Fagaceae', Cornus = 'Cornaceae', Cupressus = 'Cupressaceae', Elaeagnus = 'Elaeagnaceae', Eucalyptus = 'Myrtaceae', Fraxinus = 'Oleaceae', Juglans = 'Juglandaceae', Juniperus = 'Cupressaceae', Larix = 'Pinaceae', Liquidambar = 'Hamamelidaceae', Lithocarpus = 'Fagaceae', Malus = 'Rosaceae', Olneya = 'Fabaceae', Picea = 'Pinaceae', Pinus = 'Pinaceae', Platanus = 'Platanaceae', Populus = 'Salicaceae', Prosopis = 'Fabaceae', Prunus = 'Rosaceae', Pseudotsuga = 'Pinaceae', Quercus = 'Fagaceae', Robinia = 'Fabaceae', Salix = 'Salicaceae', Sequoia = 'Cupressaceae', Sequoiadendron = 'Cupressaceae', Taxus = 'Taxaceae', Thuja = 'Cupressaceae', Torreya = 'Taxaceae', Tsuga = 'Pinaceae', Umbellularia = 'Lauraceae')

pnw_species$family <- family_table[gsub(' ', '', pnw_species$Genus)]

# Create slash delimited file


pnw_sppout <- data.frame(Family = tolower(pnw_species$family), Genus = gsub(' ', '', pnw_species$Genus), Genus_species = gsub(' ','.',pnw_species$sciname), stringsAsFactors = F)
write.table(pnw_sppout, sep = '/', row.names = FALSE, quote = FALSE, file = file.path(fp, 'fiaspp.txt'))

# Generated tree with phylomatic using phylogeny from Zanne et al. 2014
fiatree <- read.newick(file.path(fp, 'fiatree_zanne.nwk'))
fiatree <- collapse.singles(fiatree)
fiatree$tip.label <- gsub('\\(spp\\.','\\(spp\\.\\)', fiatree$tip.label)

# Add species to the tree
spp_not_in_tree <- pnw_sppout$Genus_species[!pnw_sppout$Genus_species %in% fiatree$tip.label]
spp_not_in_tree <- gsub('\\.', ' ', spp_not_in_tree)
fiatree$tip.label <- gsub('\\.', ' ', fiatree$tip.label)

library(phytools)
for (i in spp_not_in_tree) fiatree <- add.species.to.genus(tree = fiatree, species = i, where = 'root')



##############
# Next attempt . . . 
# Qian et al 2015, Journal of Plant Ecology phylogeny (PhytoPhylo)

phytophylo <- read.tree(file.path(fp, 'S.PhyloMaker-master/PhytoPhylo'))
phymatch3 <- pnw_species$sciname %in% phytophylo$tip.label
spp_not_in_tree <- pnw_species$sciname[!phymatch3]
spp_not_in_tree <- spp_not_in_tree[!grepl('Tree',spp_not_in_tree)]


# add missing species (slow if you do on a big tree)
for (i in spp_not_in_tree) {
  phytophylo <- add.species.to.genus(tree = phytophylo, species = i, where = 'random')
  print(i)
}

save(phytophylo, file = file.path(fp,'phytophylo_taxaadded.r'))

# pare the tree down to only the species
# this could be done before adding missing species, or after
fiaphytophylo <- drop.tip(phytophylo, tip = phytophylo$tip.label[!phytophylo$tip.label %in% pnw_species$sciname])

save(fiaphytophylo, file = file.path(fp, 'phytophylo_fia.r'))