# Phylogenetic diversity metrics on FIA trees, at plot and subplot level
# Use phylogeny generated in fiaphylogeny.r
# Set up to run on cluster since it might take a while.

#fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'
fp <- '/mnt/research/nasabio'
fiapnw <- read.csv(file.path(fp, 'data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)

library(ape)
load(file.path(fp, 'data/fia/phytophylo_fia.r'))

# Get FIA species code to scientific name lookup table, and find which ones are in the tree.
fiataxa <- read.csv(file.path(fp, 'data/fia/fia_taxon_lookuptable.csv'), stringsAsFactors = FALSE)
pnw_codes <- unique(fiapnw$SPCD)
#all(pnw_codes %in% fiataxa$FIA.Code) #yes

# Convert fiapnw into a site x species matrix, at both subplot and plot level

library(dplyr)

area_by_sp <- function(dat, sppids) {
  areas <- numeric(length(sppids))
  for (i in 1:length(sppids)) {
    areas[i] <- sum(dat$basalarea[dat$SPCD == sppids[i]])
  }
  areas
}

# Subplot level
fiasums_subplot <- fiapnw %>%
  filter(STATUSCD == 1) %>%
  mutate(ba = pi * (DIA/200)^2) %>% # basal area in m2.
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT, SUBP, MEASYEAR, SPCD) %>%
  summarize(basalarea = sum(ba),
            n = n())

sppids <- sort(unique(fiasums_subplot$SPCD))

idx <- match(sppids, fiataxa$FIA.Code)
fiataxa$sciname <- paste(gsub(' ', '', fiataxa$Genus), gsub(' ', '', fiataxa$Species), sep = '_')
fiataxa$sciname[fiataxa$sciname == 'Chrysolepis_chrysophyllavar.chrysophylla'] <- 'Chrysolepis_chrysophylla'
fiataxa$sciname[fiataxa$sciname == 'Populus_balsamiferassp.Trichocarpa'] <- 'Populus_trichocarpa'


fiasubplist <- fiasums_subplot %>% do(x = area_by_sp(., sppids))
fiasubpmat <- do.call('rbind', fiasubplist$x)

# Get the species names from the lookup table that go with the numerical codes.
#sppnames <- pnw_species$sciname[match(sppids, pnw_codes)]
dimnames(fiasubpmat)[[2]] <- fiataxa$sciname[idx]

# Get rid of the unknown species.
fiasubpmat <- fiasubpmat[, dimnames(fiasubpmat)[[2]] %in% fiaphytophylo$tip.label]


# Plot level
fiasums_plot <- fiapnw %>%
  filter(STATUSCD == 1) %>%
  mutate(ba = pi * (DIA/200)^2) %>% # basal area in m2.
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT, MEASYEAR, SPCD) %>%
  summarize(basalarea = sum(ba),
            n = n())

sppids <- sort(unique(fiasums_plot$SPCD))

idx <- match(sppids, fiataxa$FIA.Code)

fiaplotlist <- fiasums_plot %>% do(x = area_by_sp(., sppids))
fiaplotmat <- do.call('rbind', fiaplotlist$x)

# Get the species names from the lookup table that go with the numerical codes.
#sppnames <- pnw_species$sciname[match(sppids, pnw_codes)]
dimnames(fiaplotmat)[[2]] <- fiataxa$sciname[idx]

# Get rid of the unknown species.
fiaplotmat <- fiaplotmat[, dimnames(fiaplotmat)[[2]] %in% fiaphytophylo$tip.label]

save(fiaplotmat, fiasubpmat, fiaplotlist, fiasubplist, file = file.path(fp, 'data/fia/fiamatrices.RData'))

##################################
# Actual calculation of PD for plot and subplot

# cophenetic distance matrix of fia phylo
fiadist <- cophenetic(fiaphytophylo)

# Calculate phylogenetic diversity metrics
library(picante)
pd_fia_plot <- pd(fiaplotmat, fiaphytophylo, include.root = TRUE)
mpd_fia_plot <- ses.mpd(fiaplotmat, fiadist, null.model = 'independentswap', abundance.weighted = TRUE, runs = 999, iterations = 1000)
mntd_fia_plot <- ses.mntd(fiaplotmat, fiadist, null.model = 'independentswap', abundance.weighted = TRUE, runs = 999, iterations = 1000)
pd_fia_subp <- pd(fiasubpmat, fiaphytophylo, include.root = TRUE)
mpd_fia_subp <- ses.mpd(fiasubpmat, fiadist, null.model = 'independentswap', abundance.weighted = TRUE, runs = 999, iterations = 1000)
mntd_fia_subp <- ses.mntd(fiasubpmat, fiadist, null.model = 'independentswap', abundance.weighted = TRUE, runs = 999, iterations = 1000)

save(fiaplotmat, fiasubpmat, fiaplotlist, fiasubplist, pd_fia_plot, mpd_fia_plot, mntd_fia_plot, pd_fia_subp, mpd_fia_subp, mntd_fia_subp, file = file.path(fp, 'data/fia/fia_pd.r'))
