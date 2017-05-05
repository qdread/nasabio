radii <- c(1000,5000,7500,10000,20000,50000,100000)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
r <- radii[task]

load('/mnt/research/nasabio/data/fia/fiaworkspace.r')

library(dplyr)

area_by_sp <- function(dat, sppids) {
  areas <- numeric(length(sppids))
  for (i in 1:length(sppids)) {
    areas[i] <- sum(dat$basalarea[dat$SPCD == sppids[i]])
  }
  areas
}

all_mats <- list()

for(p in 1:nrow(fiaalbers)) {
  if (class(fianhb_r[[p]]) == 'data.frame') {
	if (any(fianhb_r[[p]]$dist <= r)) {
		# Subset out the data frame with the nearest neighbors
		neighbs <- subset(fianhb_r[[p]], dist <= r)
		
		# Subset out the data frame with the nearest neighbors
       plotcns <- fiaalbers[c(p, neighbs$idx), ]$PLT_CN
       dat_p <- subset(fiasums_plot, PLT_CN %in% plotcns)
       # Convert into a site x species matrix
       sppids <- sort(unique(dat_p$SPCD))
       mat_p <- dat_p %>% group_by(PLT_CN) %>% do(x = area_by_sp(., sppids))
       mat_p <- do.call('rbind', mat_p$x)
	   
	   sppnames <- fiataxa$sciname[match(sppids, fiataxa$FIA.Code)]
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

save(all_mats, file = paste0('/mnt/research/nasabio/data/fia/mats/mat_',as.character(as.integer(r)),'.r'))
