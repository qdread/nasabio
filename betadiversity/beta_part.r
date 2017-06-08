# Function to run all Baselga's betapart calculations for a given matrix. This will do all pairwise calculations. Hopefully it will be fast.

# For fd and pd, requires the tree and the trait matrix in raw form, not distance matrices, unfortunately.
# It doesn't allow NA and things, so let's stick with taxonomic diversity for now.
# Can return pairwise distance matrix as well as the multisite measurement.

beta_part <- function(m, abundance=TRUE, pairwise=FALSE) {
	
	require(betapart)
	
	# Do precalculation.
	core_presence  <- betapart.core(m)
	
	# Calculate metrics.
	beta_sorensen_presence  <- unlist(beta.multi(core_presence, index.family = 'sorensen'))
	beta_jaccard_presence <- unlist(beta.multi(core_presence, index.family = 'jaccard'))
	
	if (abundance) {
		core_abundance <- betapart.core.abund(m)
		beta_bray_abundance <- unlist(beta.multi.abund(core_abundance, index.family = 'bray'))
		beta_ruzicka_abundance <- unlist(beta.multi.abund(core_abundance, index.family = 'ruzicka'))
	}
	
	# Calculate pairwise metrics if needed.
	if (pairwise) {
		beta_sorensen_presence_pair  <- beta.pair(core_presence, index.family = 'sorensen')
		beta_jaccard_presence_pair <- beta.pair(core_presence, index.family = 'jaccard')
		if (abundance) {
			beta_bray_abundance_pair <- beta.pair.abund(core_abundance, index.family = 'bray')
			beta_ruzicka_abundance_pair <- beta.pair.abund(core_abundance, index.family = 'ruzicka')
		}
	}
	# Here we can include functional and phylogenetic if needed.
	
	# Combine and return results.
	res <- data.frame(index = rep(c('sorensen', 'jaccard'), each = 3),
					  divtype = c('bal','gra','tot'),
					  abundance = FALSE,
					  beta = c(beta_sorensen_presence, beta_jaccard_presence))
	if (abundance) {
		res_abund <- data.frame(index = rep(c('bray', 'ruzicka'), each = 3),
					  divtype = c('bal','gra','tot'),
					  abundance = TRUE,
					  beta = c(beta_bray_abundance, beta_ruzicka_abundance))
		res <- rbind(res, res_abund)
	}
	
	if (!pairwise) return(res)
	res <- list(res, beta_sorensen_presence_pair, beta_jaccard_presence_pair)
	if (abundance) res <- c(res, beta_bray_abundance_pair, beta_ruzicka_abundance_pair)
	return(res)
}

# Here is a similar function using beta.div.comp from the adespatial package (does both Baselga and Podani)

beta_baselga_podani <- function(m, abundance = TRUE) {
	
	require(adespatial)
	
	beta_podani_sorensen_presence <- beta.div.comp(m, coef = 'S', quant = FALSE)$part
	beta_podani_jaccard_presence <- beta.div.comp(m, coef = 'J', quant = FALSE)$part
	beta_baselga_sorensen_presence <- beta.div.comp(m, coef = 'BS', quant = FALSE)$part
	beta_baselga_jaccard_presence <- beta.div.comp(m, coef = 'BJ', quant = FALSE)$part
	
	if (abundance) {
		beta_podani_sorensen <- beta.div.comp(m, coef = 'S', quant = FALSE)$part
		beta_podani_jaccard <- beta.div.comp(m, coef = 'J', quant = FALSE)$part
		beta_baselga_sorensen <- beta.div.comp(m, coef = 'BS', quant = FALSE)$part
		beta_baselga_jaccard <- beta.div.comp(m, coef = 'BJ', quant = FALSE)$part
	}
	
	# Combine and return results.
	res <- data.frame(family = rep(c('podani', 'baselga'), each = 10),
					  index = rep(c('sorensen', 'jaccard'), each = 5),
					  divtype = c('total','replacement','nestedness','replacement_proportion','nestedness_proportion'),
					  abundance = FALSE,
					  beta = c(beta_podani_sorensen_presence, beta_podani_jaccard_presence, beta_baselga_sorensen_presence, beta_baselga_jaccard_presence))
	
	if (abundance) {
		res_abund <- data.frame(family = rep(c('podani', 'baselga'), each = 10),
								  index = rep(c('sorensen', 'jaccard'), each = 5),
								  divtype = c('total','replacement','nestedness','replacement_proportion','nestedness_proportion'),
								  abundance = TRUE,
								  beta = c(beta_podani_sorensen, beta_podani_jaccard, beta_baselga_sorensen, beta_baselga_jaccard))
		res <- rbind(res, res_abund)		
	}
	
	return(res)
	
}