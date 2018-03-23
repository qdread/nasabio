# Find outliers in beta diversity (functional) for BBS

fd_pair_z <- bbs_betadiv_array[,,'beta_fd_pairwise_pa_z']

range(fd_pair_z,na.rm=T)

fdcut <- cut(c(fd_pair_z), breaks = seq(-7e5, 1e5, by = 100000))
table(fdcut)

fdrowmin <- apply(fd_pair_z, 1, min, na.rm = TRUE)
lower.tri(fd_pair_z, diag=FALSE)

neighbordivfromtable2 <- function(x) {
	focalpointindex <- which(bbscov_oneyear$rteNo == x$rteNo)
	neighbordists <- spDistsN1(pts = cbind(bbscov_oneyear$lon, bbscov_oneyear$lat), pt = c(x$lon, x$lat), longlat = TRUE)
	commdat <- list()
	for (i in 1:length(radii)) {
		neighbors_incircle <- bbs_betadiv_array[neighbordists <= radii[i], focalpointindex, , drop = FALSE]
		commdat[[i]] <- c(radius = radii[i], 
						  apply(neighbors_incircle[, , prop_vars, drop = FALSE], 3, median, na.rm=TRUE),
						  apply(neighbors_incircle[, , -(prop_vars), drop = FALSE], 3, median, na.rm=TRUE))
	}
	as.data.frame(do.call('rbind', commdat))
}

bbs_beta_med <- bbscov_oneyear %>%
	rowwise %>%
	do(neighbordivfromtable2(.))
	
bbs_beta_med <- cbind(bbscov_oneyear[rep(1:nrow(bbscov_oneyear), each=length(radii)),], bbs_beta_med)

# Compare
fd2ways <- bbs_beta %>% 
	select(rteNo, radius, beta_fd_pairwise_pa_z) %>%
	left_join(bbs_beta_med %>% select(rteNo, radius, beta_fd_pairwise_pa_z), by = c('rteNo','radius')) %>%
	arrange(beta_fd_pairwise_pa_z.x)