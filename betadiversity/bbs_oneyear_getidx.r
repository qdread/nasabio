# Modified 15 May: add the focal plot as row 1 of the matrix!
# Forked 25 Sep: do all years pooled to 1 year.

radii <- c(50, 75, 100, 150, 200, 300, 400, 500) * 1000
load('/mnt/research/nasabio/data/bbs/bbsworkspace_singleyear.r')

# Prep list

for (r in radii) {

print(paste('Radius', r))

all_mats <- list()
dimnames(bbscovmat_oneyear)[[1]] <- 1:nrow(bbscovmat_oneyear)

for(p in 1:nrow(bbsmat_byroute_oneyear)) {
  if (class(bbsnhb_list_oneyear[[p]]) == 'data.frame') {
	if (any(bbsnhb_list_oneyear[[p]]$dist <= r)) {
		# Subset out the data frame with the nearest neighbors
		neighbs <- subset(bbsnhb_list_oneyear[[p]], dist <= r)
		
		rowidx <- which(bbscovmat_oneyear[,1] %in% neighbs$rteNo)
		
		dat_p <- bbsmat_byroute_oneyear[c(p, rowidx), , drop = FALSE]		# Ensures that the first row of the matrix will be the focal plot.
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
	
save(all_mats, file = paste0('/mnt/research/nasabio/data/bbs/mats/oneyearroutemat_',as.character(as.integer(r)),'.r'))

}
