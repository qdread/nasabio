# Get mean, sd, and quantiles of each column of a slice

# Extract fixed and random effect draws for each iteration
# Add the fixed effect to the random effect of each region to get the slope
# Get mean, standard deviation, and quantiles for each region
# Wrap it all in a function that takes the brmsfit object as its only argument

get_stats <- function(x) {
	estimate <- apply(x, 2, mean)
	est.error <- apply(x, 2, sd)
	q025 <- apply(x, 2, quantile, probs = 0.025)
	q975 <- apply(x, 2, quantile, probs = 0.975)
	cbind(estimate, est.error, q025, q975)
}

slope_stats <- function(fit) {
	raw_fixed <- fixef(fit, summary = FALSE)
	raw_random <- ranef(fit, summary = FALSE)[[1]]
	n_reg <- dim(raw_random)[2]
	fix_plus_random <- aperm(replicate(n_reg, raw_fixed), perm = c(1,3,2)) + raw_random
	res <- list()
	for (i in 1:(dim(fix_plus_random)[3])) {
		res[[i]] <- get_stats(fix_plus_random[,,i])
	}
	res <- data.frame(effect = 'fixed_plus_random', expand.grid(dimnames(raw_random)[[2]], dimnames(raw_random)[[3]]), do.call(rbind, res))
	setNames(res, c('effect', 'region', 'parameter', 'Estimate', 'Est.Error', 'q025', 'q975'))
}