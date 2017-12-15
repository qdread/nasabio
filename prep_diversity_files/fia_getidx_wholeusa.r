# Modified 15 May: rearrange each matrix so that the focal plot is the first row of the resulting matrix.
# Forked 22 June: do in parallel for 300 only
# Forked 24 Oct: do in parallel for everything 100 and above.
# Edited 13 Dec: change input data to entire USA and make entire thing parallel.
# Edited 14 Dec: change so that neighbor data is a matrix not a data frame (less memory)

# Get radius and slice. Total 2125 tasks. Can specify how many slices per radius.
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
radii <- c(5, 10, 20, 50, 75, 100, 150, 200, 300) * 1000
n_tasks <- c(25, 25, 25, 25, 25, 500, 500, 500, 500)

combos <- data.frame(radius = rep(radii, n_tasks),
					 slice = unlist(sapply(n_tasks, function(x) 1:x)),
					 n = rep(n_tasks, n_tasks))

r <- combos$radius[task]
slice <- combos$slice[task]
n_slices <- combos$n[task]

load('/mnt/research/nasabio/data/fia/fiaworkspace_nospatial_wholeusa.r')
plotmetadata <- read.csv('/mnt/research/nasabio/data/fia/fianocoords_wholeusa.csv', stringsAsFactors = FALSE)

rowidx <- round(seq(0, nrow(plotmetadata), length.out = n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

library(dplyr)

area_by_sp <- function(dat, sppids) {
  areas <- numeric(length(sppids))
  for (i in 1:length(sppids)) {
    areas[i] <- sum(dat$basalarea[dat$SPCD == sppids[i]])
  }
  areas
}

all_mats <- list()

for(p in rowidxmin:rowidxmax){
  if (class(fianhb_r[[p]]) == 'matrix') {
	if (any(fianhb_r[[p]][,2] <= r)) {
		# Subset out the data frame with the nearest neighbors
		neighbs <- fianhb_r[[p]][fianhb_r[[p]][,2] <= r, , drop = FALSE]
		
		# Subset out the data frame with the nearest neighbors
       plotcns <- plotmetadata$PLT_CN[c(p, neighbs[,1])]
       dat_p <- subset(fiasums_plot, PLT_CN %in% plotcns)
       # Convert into a site x species matrix
       sppids <- sort(unique(dat_p$SPCD))
       mat_p <- dat_p %>% group_by(PLT_CN) %>% do(x = area_by_sp(., sppids))
	   # Sort the output so that the focal plot will be the first row of the resulting matrix.
	   focal_idx <- which(mat_p$PLT_CN == plotmetadata$PLT_CN[p])
	   mat_p <- mat_p[c(focal_idx, (1:nrow(mat_p))[-focal_idx]), ]
		
	   mat_p <- do.call('rbind', mat_p$x)
 

	   sppnames <- fiataxa$binomial_forphylo[match(sppids, fiataxa$FIA.Code)]
	   if (inherits(mat_p, 'matrix')) {
		dimnames(mat_p) <- list(1:nrow(mat_p), sppnames)
		all_mats[[length(all_mats) + 1]] <- mat_p
	  }
	  else {
		all_mats[[length(all_mats) + 1]] <- NA
	}
	}
	else {
		all_mats[[length(all_mats) + 1]] <- NA
	}
	}
  else {
	all_mats[[length(all_mats) + 1]] <- NA
  }
}

save(all_mats, file = paste0('/mnt/research/nasabio/data/fia/mats/usamat_',as.character(as.integer(r)),'_',slice,'.r'))
