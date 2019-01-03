# Modified 15 May: rearrange each matrix so that the focal plot is the first row of the resulting matrix.
# Correction 22 Aug: use workspace 2, not old workspace
# Modified 23 Oct: change to use unfuzzed coordinates
# Modified 13 Dec: change input data to entire USA.

radii <- c(5,10,20,50,75,100,150,200,300) * 1000
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
r <- radii[task]

load('/mnt/research/nasabio/data/fia/fiaworkspace_nospatial_wholeusa.r')
plotmetadata <- read.csv('/mnt/research/nasabio/data/fia/fianocoords_wholeusa.csv', stringsAsFactors = FALSE)

library(dplyr)

area_by_sp <- function(dat, sppids) {
  areas <- numeric(length(sppids))
  for (i in 1:length(sppids)) {
    areas[i] <- sum(dat$basalarea[dat$SPCD == sppids[i]])
  }
  areas
}

all_mats <- list()

for (p in 1:length(fianhb_r)) {
  if (class(fianhb_r[[p]]) == 'data.frame') {
	if (any(fianhb_r[[p]]$dist <= r)) {
		# Subset out the data frame with the nearest neighbors
		neighbs <- subset(fianhb_r[[p]], dist <= r)
		
		# Subset out the data frame with the nearest neighbors
       plotcns <- plotmetadata[c(p, neighbs$idx), ]$PLT_CN
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
	if (p%%1000 == 0) print(p)
}

save(all_mats, file = paste0('/mnt/research/nasabio/data/fia/mats/usamat_',as.character(as.integer(r)),'.r'))
