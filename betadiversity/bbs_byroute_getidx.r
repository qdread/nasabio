radii <- c(50000, 75000, 100000, 150000, 200000)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))

r <- radii[task]

load('/mnt/research/nasabio/data/bbs/bbsworkspace_byroute.r')

# Prep list

all_mats <- list()
dimnames(bbscovmat)[[1]] <- 1:nrow(bbscovmat)

for(p in 1:nrow(fixedbbsmat_byroute)) {
  if (class(bbsnhb_r[[p]]) == 'data.frame') {
	if (any(bbsnhb_r[[p]]$dist <= r)) {
		# Subset out the data frame with the nearest neighbors
		neighbs <- subset(bbsnhb_r[[p]], dist <= r)
		year_p <- bbscovmat[p,1]
		routes_p <- c(bbscovmat[p,2], neighbs$rteNo)
		# Do sequential subsetting to save time.
		yearidx <- which(bbscovmat[,1] == year_p)
		yearmat <- bbscovmat[yearidx, , drop = FALSE]
		rowidx <- as.numeric(dimnames(yearmat[which(yearmat[,2] %in% routes_p), , drop = FALSE])[[1]])
		dat_p <- fixedbbsmat_byroute[rowidx, , drop = FALSE]
		# Get rid of empty columns
		mat_p <- dat_p[, apply(dat_p, 2, sum) > 0, drop = FALSE]
		all_mats[[length(all_mats) + 1]] <- mat_p
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
	
save(all_mats, file = paste0('/mnt/research/nasabio/data/bbs/mats/routemat_',r,'.r'))
