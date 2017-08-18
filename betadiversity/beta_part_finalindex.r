# Function to run all Baselga's betapart calculations for a given matrix. This will do all pairwise calculations. Hopefully it will be fast.
# I've decided Baselga x Sorensen x Binary is the best, so that's what we're going to do.
# Later FD and PD will be incorporated.

# For fd and pd, requires the tree and the trait matrix in raw form, not distance matrices, unfortunately.
# It doesn't allow NA and things, so let's stick with taxonomic diversity for now.
# Can return pairwise distance matrix as well as the multisite measurement.

beta_part <- function(m, abundance=TRUE, pairwise=FALSE, index_family='sorensen') {
  
  require(betapart)
  m_bin <- m
  m_bin[m_bin > 0] <- 1
  
  index_family_abund <- ifelse(index_family == 'sorensen', 'bray', 'ruzicka')
  
  # Do precalculation.
  core_presence  <- betapart.core(m_bin)
  
  # Calculate metrics.
  beta_presence  <- unlist(beta.multi(core_presence, index.family = index_family))
  
  if (abundance) {
    core_abundance <- betapart.core.abund(m)
    beta_abundance <- unlist(beta.multi.abund(core_abundance, index.family = index_family_abund))
  }
  
  # Calculate pairwise metrics if needed.
  if (pairwise) {
    beta_presence_pair  <- beta.pair(core_presence, index.family = index_family)
    if (abundance) {
      beta_abundance_pair <- beta.pair.abund(core_abundance, index.family = index_family_abund)
    }
  }
  # Here we can include functional and phylogenetic if needed.
  
  # Combine and return results.
  res <- data.frame(index = rep(index_family, 3),
                    divtype = c('replacement','nestedness','total'),
                    abundance = FALSE,
                    beta = beta_presence)
  if (abundance) {
    res_abund <- data.frame(index = rep(index_family_abund, 3),
                            divtype = c('replacement','nestedness','total'),
                            abundance = TRUE,
                            beta = beta_abundance)
    res <- rbind(res, res_abund)
  }
  
  if (!pairwise) return(res)
  res <- list(res, beta_presence_pair)
  if (abundance) {
    res[[3]] <- beta_abundance_pair
  }
  return(res)
}
