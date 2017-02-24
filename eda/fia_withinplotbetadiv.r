# Within-plot beta-diversity for FIA

fp <- '/mnt/research/nasabio'
fiapnw <- read.csv(file.path(fp, 'data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)

library(ape)
load(file.path(fp, 'data/fia/phytophylo_fia.r'))
fiadist <- cophenetic(fiaphytophylo)

fiataxa <- read.csv(file.path(fp, 'data/fia/fia_taxon_lookuptable.csv'), stringsAsFactors = FALSE)
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

pnw_species <- subset(pnw_species, !grepl('Tree', sciname)) # unidentified species
#pnw_species$sciname <- gsub('_', ' ', pnw_species$sciname)
pnw_species <- pnw_species[order(pnw_species$sciname), ]

# Calculate basal area at subplot level
library(dplyr)

fiasums_subplot <- fiapnw %>%
  filter(STATUSCD == 1) %>%
  mutate(ba = pi * (DIA/200)^2) %>% # basal area in m2.
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT, SUBP, MEASYEAR, SPCD) %>%
  summarize(basalarea = sum(ba),
            n = n())


makecommmat <- function(dat) {
  sppids <- sort(unique(dat$SPCD))
  mat <- dat %>% group_by(SUBP) %>% do(x = sapply(sppids, function(z) sum(.$basalarea[.$SPCD == z])))
  mat <- do.call('rbind', mat$x)
  if(!is.null(mat)) {
    if(nrow(mat) > 1) {
      # Fix the species names to match the phylogeny, and get rid of the unknown species.
      sppnames <- pnw_species$sciname[match(sppids, pnw_codes)]
      dimnames(mat)[[1]] <- 1:nrow(mat)
      dimnames(mat)[[2]] <- sppnames
      mat <- mat[, dimnames(mat)[[2]] %in% fiaphytophylo$tip.label, drop = FALSE]
    }
  }
  mat
}

# Within-plot beta diversity.
# Generate matrices:
fiamats <- fiasums_subplot %>% 
  ungroup %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT, MEASYEAR) %>%
  do(comm = makecommmat(.))

