radii <- c(1000,2000,3000,4000,5000)
task <- as.numeric(Sys.getenv('PBS_ARRAYID'))
n_slices <- 50

combos <- data.frame(slice = rep(1:n_slices, times = length(radii)), radius = rep(radii, each = n_slices))

r <- combos$radius[task]
slice <- combos$slice[task]


load('/mnt/research/nasabio/data/bbs/bbsworkspace.r')

# Prep list

all_mats <- list()
dimnames(bbscovmat)[[1]] <- 1:nrow(bbscovmat)

# Determine row indices for the slice of the matrix to be used.
rowidx <- round(seq(0,nrow(fixedbbsmat),length.out=n_slices + 1))
rowidxmin <- rowidx[slice]+1
rowidxmax <- rowidx[slice+1]

for(p in rowidxmin:rowidxmax) {
  if (class(bbsnhb_r[[p]]) == 'data.frame') {
	if (any(bbsnhb_r[[p]]$dist <= r)) {
		# Subset out the data frame with the nearest neighbors
		neighbs <- subset(bbsnhb_r[[p]], dist <= r)
		year_p <- bbscovmat[p,1]
		route_p <- bbscovmat[p,2]
		stops_p <- c(bbscovmat[p,3], neighbs$Stop)
		# Do sequential subsetting to save time.
		yearidx <- which(bbscovmat[,1] == year_p)
		yearmat <- bbscovmat[yearidx, , drop = FALSE]
		routeidx <- which(yearmat[,2] == route_p)
		routemat <- yearmat[routeidx, , drop = FALSE]
		rowidx <- as.numeric(dimnames(routemat[which(routemat[,3] %in% stops_p), , drop = FALSE])[[1]])
		dat_p <- fixedbbsmat[rowidx, , drop = FALSE]
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
	
save(all_mats, file = paste0('/mnt/research/nasabio/data/bbs/mats/mat_',r,'_',slice,'.r'))
