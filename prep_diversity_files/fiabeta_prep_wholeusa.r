# FIA taxonomic and phylogenetic beta diversity calculation for the cluster.
# Forked 17 April. Just preparatory work, actual diversity calculated in a different script.
# Modified 17 May: add rows to matrix that have zero species.
# Modified 22 Aug: add better imputed traits
# Modified 20 Oct: new coordinates (workspace should now only contain pointers to the true coordinates)
# New forked version created 13 Dec: do the entire USA.
# Modified 14 Dec: use spDistsN1() which is faster than dnearneigh()
# Modified 26 Nov 2018: get rid of anything with TPA_UNADJ identifying it as not being in a subplot.
# Modified 27 Nov 2018: instead of using TPA_UNADJ to identify trees to remove, use the CSV with distance to plot centers. Also remove the calculation of neighbors here since it is done elsewhere.
# Modified 12 Dec 2018: add plantation status to workspace

fp <- '/mnt/research/nasabio'
fiaall <- read.csv(file.path(fp, 'data/fia/treedata10nov/finley_trees_continental_US_most_recent_evaluations_nov8_2017.csv'), stringsAsFactors = FALSE)

#phylogenetic distance matrix
library(ape)
load(file.path(fp, 'data/fia/treedata10nov/phylogenies_allfia.r'))
fiadist <- cophenetic(fullphylo)

# Get FIA species code to scientific name lookup table
fiataxa <- read.csv(file.path(fp, 'data/fia/treedata10nov/lookup_table_allfia.csv'), stringsAsFactors = FALSE)
allfia_codes <- unique(fiaall$SPCD)

# Modification 26 Nov 2018: here, get rid of the macroplots by filtering on TPA_UNADJ (microplots are okay)
fiadists <- read.csv(file.path(fp, 'data/fia/finley_pnw_tree_distances_to_subp_center_nov26_2018.csv'))
trees_in_macroplots <- with(fiadists, TREE_CN[DIST > 24])
fiaall <- fiaall[!fiaall$TREE_CN %in% trees_in_macroplots, ]

# For all the FIA taxa that are represented by two species codes after trait and phylo QC have been run,
# assign the single correct species code to them.
# Write this to work in the general case instead of as crappy ad hoc code.

id_table <- table(fiataxa$binomial_forphylo)
duplicate_code_rows <- which(fiataxa$binomial_forphylo %in% names(id_table[id_table > 1]))
good_names <- unique(fiataxa$binomial_forphylo[duplicate_code_rows])
bad_names <- fiataxa$binomial[duplicate_code_rows][!fiataxa$binomial[duplicate_code_rows] %in% good_names]
good_names_ordered <- fiataxa$binomial_forphylo[duplicate_code_rows][!fiataxa$binomial[duplicate_code_rows] %in% good_names]
good_codes <- fiataxa$FIA.Code[match(good_names_ordered, fiataxa$binomial)]
bad_codes <- fiataxa$FIA.Code[match(bad_names, fiataxa$binomial)]
names(good_codes) <- good_names_ordered # just to confirm they line up.
names(bad_codes) <- bad_names

# Replace the bad codes with good codes (merging invalid species with the valid one).
for (i in 1:length(good_codes)) {
	fiaall$SPCD[fiaall$SPCD == bad_codes[i]] <- good_codes[i]
}

# Convert fiaall into a site x species matrix at plot level

library(dplyr)

area_by_sp <- function(dat, sppids) {
  areas <- numeric(length(sppids))
  for (i in 1:length(sppids)) {
    areas[i] <- sum(dat$basalarea[dat$SPCD == sppids[i]])
  }
  areas
}

# Plot level
fiasums_plot <- fiaall %>%
  filter(STATUSCD == 1) %>%
  mutate(ba = pi * (DIA/200)^2) %>% # basal area in m2.
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT, MEASYEAR, SPCD) %>%
  summarize(basalarea = sum(ba),
            n = n())
			
sppids <- sort(unique(fiasums_plot$SPCD))

idx <- match(sppids, fiataxa$FIA.Code)

fiacoords <- read.csv('/mnt/home/qdr/data/allfia.csv', stringsAsFactors = FALSE)
fiacoords <- fiacoords %>%
	filter(!is.na(ACTUAL_LAT), CN %in% fiaall$PLT_CN) %>%
	rename(PLT_CN = CN, lat = ACTUAL_LAT, lon = ACTUAL_LON)

# The summing and reshaping of area by species takes a lot longer with the full dataset.
fiaplotlist <- list()
pb <- txtProgressBar(0, nrow(fiacoords), style=3)
for (i in 1:nrow(fiacoords)) {
	fiaplotlist[[i]] <- area_by_sp(subset(fiasums_plot, PLT_CN == fiacoords$PLT_CN[i]), sppids)
	setTxtProgressBar(pb,i)
}
close(pb)
fiaplotmat <- do.call('rbind', fiaplotlist)

# Get the species names from the lookup table that go with the numerical codes.
dimnames(fiaplotmat)[[2]] <- fiataxa$binomial_forphylo[idx]

# Get rid of the unknown species.
fiaplotmat <- fiaplotmat[, dimnames(fiaplotmat)[[2]] %in% fullphylo$tip.label]


# Reproject FIA plots and calculate the distance in m from each plot to all other plots.
library(rgdal)

aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'

fiaspatial <- SpatialPointsDataFrame(coords = data.frame(x=fiacoords$lon, y=fiacoords$lat),
                                     data = data.frame(PLT_CN = fiacoords$PLT_CN),
                                     proj4string = CRS(wgs_crs)
)

fiaalbers <- spTransform(fiaspatial, CRSobj = CRS(aea_crs))

# Generate distance matrix for functional beta-diversity
traits_imputed <- read.csv(file.path(fp,'data/fia/treedata10nov/traits_imputed_allfia.csv'), stringsAsFactors = FALSE, row.names = 1)
library(FD)

trydist <- gowdis(scale(traits_imputed))

# Edit 29 Nov 2018: create a table with the state codes for lookup
fia_statecodes <- fiacoords[, 'PLT_CN', drop = FALSE] %>% left_join(unique(fiaall[, c('PLT_CN', 'STATECD', 'COUNTYCD')]))

# Edit 12 Dec 2018: add plantation status
plantation <- read.csv('/mnt/research/nasabio/data/fia/plotcond/plantation.csv')
fiaplantation <- left_join(fiacoords, plantation) %>% dplyr::select(PLT_CN, plantation)

# New saved workspace that does not actually contain the coordinates.
save(fiadist, trydist, traits_imputed, fullphylo, fiataxa, fiasums_plot, sppids, fiaplotmat, fiaplantation, file = '/mnt/research/nasabio/data/fia/fiaworkspace_nospatial_wholeusa_2018.r')
save(fiaspatial, fiaalbers, fiacoords, file = '/mnt/home/qdr/data/fiaworkspace_spatial_wholeusa_2018.r')
write.csv(fiaalbers@data, file = '/mnt/research/nasabio/data/fia/fianocoords_wholeusa_2018.csv', row.names = FALSE)
write.csv(fia_statecodes, file = '/mnt/research/nasabio/data/fia/fiastatecodes.csv', row.names = FALSE)