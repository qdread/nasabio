# Test beta_part code with functional diversity and phylogenetic diversity
# Also try to deal with issues where a single faulty plot will make the multisite index NA.

library(betapart)

# Load some FIA data to test
load('/mnt/research/nasabio/data/fia/mats/newmat_100000.r')
load('/mnt/research/nasabio/data/fia/fiaworkspace2.r') 


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
    core_func <- functional.betapart.core(m_bin, traits = trait_mat)
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

######################################
# Test

set.seed(444)
idx_test <- sample(length(all_mats), 10)

fia_list <- list()

pb <- txtProgressBar(0, length(idx_test), style=3)

for (p in 1:length(idx_test)) {
  setTxtProgressBar(pb, p)
  mat_p <- all_mats[[idx_test[p]]]
  
  if(inherits(mat_p, 'matrix')) {
    if(nrow(mat_p) > 1 & ncol(mat_p) > 1) {
      
      # Fix the species names to match the phylogeny, and get rid of the unknown species.
      
      mat_p <- mat_p[, dimnames(mat_p)[[2]] %in% pnwphylo$tip.label, drop = FALSE]
      beta_td <- beta_part(m = mat_p, abundance = TRUE, pairwise = FALSE, index_family = 'sorensen')
      
      fia_list[[p]] <- cbind(nneighb = nrow(mat_p) - 1, beta_td)
      
    }
    else {
      fia_list[[p]] <- cbind(nneighb = NA, null_result)
    }
  }
  else {
    fia_list[[p]] <- cbind(nneighb = NA, null_result)
  }
}

close(pb)


##############################################
# try with other diversity types.

pb <- txtProgressBar(0, length(idx_test), style=3)

for (p in 1:length(idx_test)) {
  setTxtProgressBar(pb, p)
  mat_p <- all_mats[[idx_test[p]]]
  
  if(inherits(mat_p, 'matrix')) {
    if(nrow(mat_p) > 1 & ncol(mat_p) > 1) {
      
      # Fix the species names to match the phylogeny, and get rid of the unknown species.
      
      mat_p <- mat_p[, dimnames(mat_p)[[2]] %in% pnwphylo$tip.label, drop = FALSE]
      beta_td <- beta_part(m = mat_p, abundance = TRUE, pairwise = FALSE, index_family = 'sorensen', FD=TRUE, PD=TRUE, phylo_tree=pnwphylo)

      fia_list[[p]] <- cbind(nneighb = nrow(mat_p) - 1, beta_td)
      
    }
    else {
      fia_list[[p]] <- cbind(nneighb = NA, null_result)
    }
  }
  else {
    fia_list[[p]] <- cbind(nneighb = NA, null_result)
  }
}

close(pb)