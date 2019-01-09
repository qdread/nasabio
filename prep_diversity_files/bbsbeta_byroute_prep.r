# Beta-diversity calculation for breeding bird survey data
# Run in parallel on HPCC (array job with a different radius for each task)
# QDR, 07 Mar 2017
# Project: NASABioXGeo

# Modification 13 Dec 2018: Go back to raw data! Drastically simplify how the species are assigned. Just plain get rid of the hybrids, unknowns, and assign subspecies to "parent" species.
# Modification 19 Apr 2018: Replace centroid coordinates with midpoints of the segments (gets rid of some routes but is more accurate)
# Modification 22 Sep 2017: (1) phylo distance matrix based on a single consensus tree of 1000 samples, (2) make single presence-absence for entire time period 2007-2016.
# Updated on 31 May 2017: new 2016 routes+coordinates.
# Modification 05 May 2017: New coordinates (correct)
# Rereforked on 18 April 2017: Do by route.
# Modification 11 April 2017: taking way too long so create an array job by radius.
# Reforked on 10 April 2017: Only preparatory work. Save "chunks" to send small pieces to HPCC.
# Forked on 06 April 2017: do taxonomic diversity only, which makes it quicker. We can just do 1 task per radius now.
# Modification 08 Mar 2017: Make doubly parallel by giving a different radius to each task
# but then also splitting the loop into 10 chunks since it will take forever to run otherwise.

# 1. Define community either as aggregated route or by stop

# (Both are done, for flavors MS it's aggregated by route, for scaling MS it's by stop)

# 2. Create community matrix.

library(dplyr)
library(reshape2)

bbsspp <- read.csv('/mnt/research/nasabio/data/bbs/bbs_species_lookup_table_modified.csv', stringsAsFactors = FALSE)

fp <- '/mnt/research/aquaxterra/DATA/raw_data/BBS/DataFiles26may2017/50-StopData/1997ToPresent_SurveyWide'

bbsdf <- list()

for (i in 1:10) {
	bbsdf[[i]] <- read.csv(file.path(fp, paste0('fifty',i,'.csv')))
}

bbsdf <- do.call('rbind', bbsdf)

# Use only surveys done with the normal protocol.
bbsdf <- subset(bbsdf, RPID == 101)

# Combine state number and route number into single rteNo
rteNo <- character(nrow(bbsdf))
for (i in 1:nrow(bbsdf)) {
	rteNo[i] <- with(bbsdf, paste(statenum[i], paste(rep('0', 3-nchar(as.character(Route[i]))), collapse=''), Route[i], sep = ''))
}

bbsdf$rteNo <- rteNo

# Convert to long format
bbslong <- melt(bbsdf, id.vars = grep('Stop', names(bbsdf), invert=TRUE), variable.name = 'Stop', value.name = 'n')
bbslong <- bbslong %>% filter(n>0) # Get rid of zero rows.

# Use species lookup table to consolidate the bbs matrix, getting rid of hybrids and unknowns and including subspecies with their parent species
bbs_aous <- unique(bbslong$AOU)
table(bbs_aous %in% bbsspp$AOU)

# let's do the following:
# (1) Anything that is a hybrid, get rid of. (7) This totals a vanishingly tiny percent of the total. (0.0002%)
# (2) Anything that is a subspecies assign to parent species (15)
# (3) Anything that is an unknown species within a genus, or an unknown genus, get rid of. They total 0.1% of individuals.

unkaou <- bbsspp$AOU[bbsspp$Type %in% c('unknown_genus', 'unknown_sp')]
unkaou <- c(unkaou, 6960) # Removes another unknown species only present in Alaska
unktotal <- sum(bbslong$n[bbslong$AOU %in% unkaou])
unktotal/sum(bbslong$n)
hybridaou <- bbsspp$AOU[bbsspp$Type %in% c('hybrid')]
sum(bbslong$n[bbslong$AOU %in% hybridaou])/sum(bbslong$n)

# Get rid of the ones that need to be gotten rid of.
aous_exclude <- c(unkaou, hybridaou)
bbslong_cleaned <- bbslong %>% filter(!AOU %in% aous_exclude)

# Create a mapping of species to their parent taxa. 
sspaou <- bbsspp$AOU[bbsspp$Type %in% 'subspecies']
parentaou <- as.numeric(bbsspp$AOU_list[bbsspp$Type %in% 'subspecies'])

# Must include the scrub jays still as a manual correction.
sspaou <- c(sspaou, 4812, 4813)
parentaou <- c(parentaou, 4811, 4811)

ssptotal <- sum(bbslong$n[bbslong$AOU %in% sspaou])
ssptotal/sum(bbslong$n)


# Replace subspecies AOUs with the one from their parent species
bbslong_cleaned$AOU[bbslong_cleaned$AOU %in% sspaou] <- na.omit(parentaou[match(bbslong_cleaned$AOU, sspaou)])

# Do final corrections mapping the species not present in the phylogeny to parent taxa they were recently split from.
recentlysplitaou <- c(4172, 7222, 16600, 5738)
recentlysplitparentaou <- c(4171, 7221, 7180, 5739)
recentlysplittotal <- sum(bbslong$n[bbslong$AOU %in% recentlysplitaou])
recentlysplittotal/sum(bbslong$n)

# Replace "recently split" AOUs with the one from their related "parent" species
bbslong_cleaned$AOU[bbslong_cleaned$AOU %in% recentlysplitaou] <- na.omit(recentlysplitparentaou[match(bbslong_cleaned$AOU, recentlysplitaou)])

sppids <- sort(unique(bbslong_cleaned$AOU)) # 631 final species.

# Convert to a site by species matrix (site is a route by year combination)
get_spp <- function(dat, sppids) {
  ns <- dat$n
  idx <- match(dat$AOU, sppids)
  res <- rep(0, length(sppids))
  res[idx] <- ns
  res
}

bbsmat_byroute <- bbslong_cleaned %>% group_by(year, rteNo) %>% do(s = get_spp(., sppids=sppids))
bbsgrps_byroute <- select(bbsmat_byroute, year, rteNo)
bbsmat_byroute <- do.call('rbind', bbsmat_byroute$s)

# 2b. Create community matrix by stop as well (added 09 Jan 2019)

bbsmat_bystop <- bbslong_cleaned %>% group_by(year, rteNo, Stop) %>% do(s = get_spp(., sppids=sppids)) # Takes 10 minutes on HPCC.
bbsgrps_bystop <- select(bbsmat_bystop, year, rteNo, Stop)
bbsmat_bystop <- do.call('rbind', bbsmat_bystop$s)
save(bbsgrps_bystop, bbsmat_bystop, file = '/mnt/research/nasabio/data/bbs/bbsmat_bystop.RData')

# 3. Project lat-long to Albers.

#bbsalbers <- read.csv('/mnt/research/nasabio/data/bbs/bbs_correct_route_centroids.csv') # Old centroids. Can still use if desired.
bbsalbers <- read.csv('/mnt/research/nasabio/data/bbs/bbs_route_midpoints.csv')

# 4. Create cophenetic distance matrix from bird phylogeny. Use consensus tree calculatd from 1000 trees randomly sampled from the posterior distribution of 10000 trees.
# (Update 22 sep 2017)

library(ape)
eric_cons_tree <- read.tree('/mnt/research/nasabio/data/bbs/ericson_cons.tre')
ericdist <- cophenetic(eric_cons_tree)

# 5. Create trait distance matrix (Gower) from bird traits.

birdtrait <- read.csv('/mnt/research/aquaxterra/DATA/raw_data/BBS/bird_traits/birdtraitmerged.csv', stringsAsFactors = FALSE)

library(FD)

birdtrait[birdtrait == -999] <- NA
nocturnalbirdAOUs <- birdtrait$AOU[birdtrait$Nocturnal == 1]
birdtrait_diurnal <- subset(birdtrait, Nocturnal != 1 | is.na(Nocturnal))
traitnames <- names(birdtrait)[c(15:24, 29:36, 46:50, 53, 55:59)]

birdtraitdist <- as.matrix(gowdis(birdtrait_diurnal[, traitnames]))
dimnames(birdtraitdist)[[1]] <- dimnames(birdtraitdist)[[2]] <- birdtrait_diurnal$Latin_Name

# 5b. Perform any last data cleaning that is necessary.

# Remove stops that don't have coordinates from the dataset.
# Here check why some of the coordinates are missing.

bbsgrps_byroute <- bbsgrps_byroute %>% mutate(rteNo=as.numeric(rteNo)) %>% left_join(bbsalbers, by = 'rteNo')
has_coords <- !is.na(bbsgrps_byroute$lon)

# Set dimnames of community matrix, trait distance matrix, and phylogenetic distance matrix to match.

dimnames_matrix <- rep(NA, length(sppids))
for (i in 1:length(sppids)) {
	names_i <- bbsspp[bbsspp$AOU == sppids[i], c('Latin_Name')]
	if (length(names_i)>0) dimnames_matrix[i] <- names_i[1]
}

dimnames(bbsmat_byroute)[[2]] <- dimnames_matrix

tlabel1 <- eric_cons_tree$tip.label
tlabel1 <- gsub('_', ' ', tlabel1)

phymatchidx <- rep(NA,length(tlabel1))

for (i in 1:length(tlabel1)) {
	phymatchidx[i] <- c(which(bbsspp$Latin_Name_clean == tlabel1[i] | bbsspp$Latin_Name_synonym == tlabel1[i] | bbsspp$Latin_Name_synonym2 == tlabel1[i]), NA)[1]
}
tlabelaou <- bbsspp$AOU[phymatchidx]

dimnames_tlabel <- rep(NA, length(tlabelaou))
for (i in 1:length(tlabelaou)) {
	names_i <- bbsspp[bbsspp$AOU == tlabelaou[i], c('Latin_Name')]
	if (length(names_i)>0) dimnames_tlabel[i] <- names_i[1]
}

dimnames(ericdist)[[1]] <- dimnames(ericdist)[[2]] <- dimnames_tlabel

# Remove nocturnal species, rows without coordinates, and rows and columns with zero sum.
ns <- colSums(bbsmat_byroute)
rs <- rowSums(bbsmat_byroute)
nocturnalbirds <- birdtrait$Latin_Name[birdtrait$Nocturnal == 1]
fixedbbsmat_byroute <- bbsmat_byroute[has_coords & rs != 0, !(dimnames(bbsmat_byroute)[[2]] %in% nocturnalbirds) & ns != 0]

#bbsalbers <- bbsalbers[has_coords & rs != 0, ]
bbsgrps_byroute <- bbsgrps_byroute[has_coords & rs != 0, ]

##################################
# Section to keep all the rows, even the ones without midpoints.
fixedbbsmat_byroute_allrows <- bbsmat_byroute[rs != 0, !(dimnames(bbsmat_byroute)[[2]] %in% nocturnalbirds) & ns != 0]

#bbsalbers <- bbsalbers[has_coords & rs != 0, ]
bbsgrps_byroute_allrows <- bbsgrps_byroute[rs != 0, ]

# Put together as a matrix.
bbsmat_allrows <- as.matrix(cbind(bbsgrps_byroute_allrows[,1:2], fixedbbsmat_byroute_allrows))
save(bbsmat_allrows, file = '~/bbsmatrix_byroute.RData')
##################################

# Subset the distance matrices by the species that occur in the BBS data

final_spnames <- dimnames(fixedbbsmat_byroute)[[2]]
ericdist <- ericdist[dimnames(ericdist)[[1]] %in% final_spnames, dimnames(ericdist)[[2]] %in% final_spnames]
birdtraitdist <- birdtraitdist[dimnames(birdtraitdist)[[1]] %in% final_spnames, dimnames(birdtraitdist)[[2]] %in% final_spnames]

save(ericdist, birdtraitdist, file = '/mnt/research/nasabio/data/bbs/bbspdfddist.r')


# 6. For the given radius, get pairwise distances among stops and the list of all neighbors.

# Fast function to get neighbor distances
# Refer to http://gis.stackexchange.com/questions/132384/distance-to-nearest-point-for-every-point-same-spatialpointsdataframe-in-r
getNeighbors <- function(dat, radius) {
  library(spdep)
  coords <- cbind(dat$lon_aea, dat$lat_aea)
  idlist <- dnearneigh(coords, 0, radius)
  distlist <- nbdists(idlist, coords)
  dflist <- list()
  for (i in 1:length(idlist)) {
    if (any(distlist[[i]] <= radius)) {
      dflist[[i]] <- data.frame(idx = idlist[[i]], rteNo = dat$rteNo[idlist[[i]]], dist = distlist[[i]][distlist[[i]] <= radius])
    }
    else {
      dflist[[i]] <- NA
    }
  }
  dflist
}

# Calculate distances and identities of all neighbors within the maximum radius
names(bbsgrps_byroute) <- c('year','rteNo','lon','lat','lon_aea','lat_aea')
bbscov <- bbsgrps_byroute

# For optimization purposes, convert covariates to a matrix.
bbscovmat <- as.matrix(bbscov)


bbsnhb_list <- as.data.frame(bbscovmat) %>% group_by(year) %>% do(l = getNeighbors(., radius = 5e5)) 
# Flatten this into one list
bbsnhb_r <- do.call('c', bbsnhb_list$l)

save(bbsnhb_r, bbscov, bbscovmat, fixedbbsmat_byroute, file = '/mnt/research/nasabio/data/bbs/bbsworkspace_byroute.r')
write.csv(fixedbbsmat_byroute, file = '/mnt/research/nasabio/data/bbs/bbs_plot_matrix.csv', row.names = FALSE)

######
# Combine 2007-2016 into a single year.

consolidate_years <- function(x) {
	mat_x <- fixedbbsmat_byroute[x$rowidx, , drop = FALSE]
	as.numeric(apply(mat_x, 2, sum) > 0)
}

bbs_consol <- bbscov %>%
	mutate(rowidx = 1:nrow(bbscov)) %>%
	filter(year >= 2007 & year <= 2016) %>%
	group_by(rteNo, lon, lat, lon_aea, lat_aea) %>%
	do(x = consolidate_years(.))
	
bbsmat_byroute_oneyear <- do.call('rbind', bbs_consol$x)
dimnames(bbsmat_byroute_oneyear)[[2]] <- dimnames(fixedbbsmat_byroute)[[2]]

bbscov_oneyear <- bbs_consol %>% select(-x)
bbscovmat_oneyear <- as.matrix(bbscov_oneyear)

bbsnhb_list_oneyear <- getNeighbors(dat = as.data.frame(bbscovmat_oneyear), radius = 5e5)

save(bbsnhb_list_oneyear, bbscov_oneyear, bbscovmat_oneyear, bbsmat_byroute_oneyear, file = '/mnt/research/nasabio/data/bbs/bbsworkspace_singleyear.r')