# Function to run all Baselga's betapart calculations for a given matrix. This will do all pairwise calculations. Hopefully it will be fast.
# I've decided Baselga x Sorensen x Binary is the best, so that's what we're going to do.
# TD is done now along with PD (binary only). Later FD will be incorporated.

# For fd and pd, requires the tree and the trait matrix in raw form, not distance matrices, unfortunately.
# PD works, but FD has too many limitations. It requires the number of species at each site to always exceed the number of traits used.
# That is too restrictive so we will have to figure out a better way of doing FD later.
# Can return pairwise distance matrix as well as the multisite measurement.

beta_part <- function(m, abundance=TRUE, pairwise=FALSE, index_family='sorensen', TD=TRUE, FD=FALSE, PD=FALSE, trait_mat = NULL, phylo_tree = NULL) {
  
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
  
  if (FD) {
    trait_mat <- trait_mat[dimnames(trait_mat)[[1]] %in% dimnames(m_bin)[[2]], ]
    core_func <- functional.betapart.core(m_bin, traits = as.matrix(trait_mat))
    beta_func <- unlist(functional.beta.multi(core_func, index.family = index_family))
  }
  
  if (PD) {
    core_phy <- phylo.betapart.core(m_bin, tree = phylo_tree)
    beta_phy <- unlist(phylo.beta.multi(core_phy, index.family = index_family))
  }
  
  # Calculate pairwise metrics if needed.
  if (pairwise) {
    beta_presence_pair  <- beta.pair(core_presence, index.family = index_family)
    if (abundance) {
      beta_abundance_pair <- beta.pair.abund(core_abundance, index.family = index_family_abund)
    }
    if (FD) {
      beta_func_pair <- functional.beta.pair(core_func, traits = trait_mat, index.family = index_family)
    }
    if (PD) {
      beta_phylo_pair <- phylo.beta.pair(core_phy, tree = phylo_tree, index.family = index_family)
    }
  }

  # Combine and return results.
  res <- data.frame(index = rep(index_family, 3),
                    diversity = 'taxonomic',
                    partition = c('replacement','nestedness','total'),
                    abundance = FALSE,
                    beta = beta_presence)
  if (abundance) {
    res_abund <- data.frame(index = rep(index_family_abund, 3),
                            diversity = 'taxonomic',
                            partition = c('replacement','nestedness','total'),
                            abundance = TRUE,
                            beta = beta_abundance)
    res <- rbind(res, res_abund)
  }
  if (FD) {
    res_func <- data.frame(index = rep(index_family, 3),
                           diversity = 'functional',
                           partition = c('replacement','nestedness','total'),
                           abundance = FALSE,
                           beta = beta_func)
    res <- rbind(res, res_func)
  }
  if (PD) {
    res_phy <- data.frame(index = rep(index_family, 3),
                          diversity = 'phylogenetic',
                          partition = c('replacement','nestedness','total'),
                          abundance = FALSE,
                          beta = beta_phy)
    res <- rbind(res, res_phy)
  }
  
  if (!pairwise) return(res)
  res <- list(res, beta_presence_pair)
  if (abundance) res[[length(res) + 1]] <- beta_abundance_pair
  if (FD) res[[length(res) + 1]] <- beta_func_pair
  if (PD) res[[length(res) + 1]] <- beta_phy_pair
  return(res)
}