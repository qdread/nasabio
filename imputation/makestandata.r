# Z is design matrix with dimensions (n trait * n indiv) x (n trait * n species) to make sure each individual maps to the right species

# Define function to make Z

make_Z <- function(traitdat) {
  Z <- matrix(0, nrow = nrow(traitdat), ncol = length(unique(traitdat$trait)) * length(unique(traitdat$species)))
  
  # For each row in traitdat, put a 1 in the correct row and column of Z.
  species_trait_combo <- with(traitdat, paste(species,trait))
  col_idx <- match(species_trait_combo, unique(species_trait_combo))
  for (i in 1:length(col_idx)) Z[i, col_idx[i]] <- 1
  return(Z)
}

# X is predictor matrix with dimensions (n trait * n indiv) x (n trait * n predictors)

# Define function to make X (Jay is writing)

make_X <- function(preds, m) {
  preds <- as.matrix(preds)
  p = ncol(preds) 
  n = nrow(preds)
  #X <- matrix(nrow=(m*nrow(preds)),ncol = m*(p+1), byrow = TRUE)
  bigList <- list()
  for(k in 1:n){
    mlist <- list()
    for(i in 1:m){
      
      X <- matrix(0,nrow=m,ncol=p+1) 
      X[i,] <- c(1,as.numeric(preds[k,]))
      mlist[[length(mlist) + 1]] <- X
    }
    ilist <- do.call(cbind,mlist)
    bigList[[k]] <- ilist
  }
  return(do.call(rbind,bigList))
}

# Y is vector of trait values by individuals/species, all strung into one long vector

# Define function to make Y

make_Y <- function(traitdat) {
  # traitdata has first column as species, rest as traits
  Y <- do.call(c, traitdat[,-1])
  names(Y) <- paste(traitdat$species, rep(names(traitdat)[-1], each = nrow(traitdat)), sep = '_')
  return(Y)
}

# Define function to make data list for stan model

make_standatalist <- function(traitdat, predictors, phy, evolution_model = c('brownian','ou')) {
  require(reshape2)
  require(Matrix)
  traitdat_melted <- melt(traitdat, id.vars = 'species', variable.name = 'trait')
  
  splist <- unique(traitdat$species)
  
  # R: phylogenetic VCV matrix. Can be either Brownian or Ornstein.
  if (evolution_model[1] == 'brownian') {
    R <- as.matrix(vcv(phy))
  }
  if (evolution_model[1] == 'ou') {
    # Function for vcv under OU model is in PIGShift
    library(PIGShift)
    R <- OU.vcv(phy, theta = 1) ### Need to check whether this is sensitive to our choice of theta.
    dimnames(R) <- list(phy$tip.label, phy$tip.label) # Label matrix rows and columns with species names.
  }
  
  R <- R[splist, splist] # Sort matrix so that it's in the same order as the trait and environment data frames.
  
  # Force R to be symmetric
  R <- as.matrix(forceSymmetric(R, "U"))
  
  m <- length(unique(traitdat_melted$trait))
  X <- make_X(predictors[,-1], m)
  Z <- make_Z(traitdat_melted)
  Y <- make_Y(traitdat)
  
  datlist <- list(m = m,
                  n = length(unique(traitdat$species)),
                  N = nrow(traitdat_melted)/m,
                  p = ncol(predictors), # this only works out because predictors has 1 extra column for species, then we add 1 to p for the intercept
                  X = X,
                  Y = Y,
                  Z = Z,
                  R = R)
  
  return(datlist)
}

make_standatalist_missing <- function(traitdat, predictors, phy, evolution_model = c('brownian','ou')) {
  require(reshape2)
  require(Matrix)
  traitdat_melted <- melt(traitdat, id.vars = 'species', variable.name = 'trait')
  
  splist <- unique(traitdat$species)
  
  # R: phylogenetic VCV matrix. Can be either Brownian or Ornstein.
  if (evolution_model[1] == 'brownian') {
    R <- as.matrix(vcv(phy))
  }
  if (evolution_model[1] == 'ou') {
    # Function for vcv under OU model is in PIGShift
    library(PIGShift)
    R <- OU.vcv(phy, theta = 1) ### Need to check whether this is sensitive to our choice of theta.
    dimnames(R) <- list(phy$tip.label, phy$tip.label) # Label matrix rows and columns with species names.
  }
  
  R <- R[splist, splist] # Sort matrix so that it's in the same order as the trait and environment data frames.
  
  # Force R to be symmetric
  R <- as.matrix(forceSymmetric(R, "U"))
  
  m <- length(unique(traitdat_melted$trait))
  X <- make_X(predictors[,-1], m)
  Z <- make_Z(traitdat_melted)
  Y <- make_Y(traitdat)
  
  index_obs <- which(!is.na(Y))
  index_mis <- which(is.na(Y))
  Y_obs <- Y[index_obs]
  N_obs <- length(index_obs)
  N_mis <- length(index_mis)
  
  datlist <- list(m = m,
                  n = length(unique(traitdat$species)),
                  N = nrow(traitdat_melted)/m,
                  p = ncol(predictors), # this only works out because predictors has 1 extra column for species, then we add 1 to p for the intercept
                  X = X,
                  Y_obs = Y_obs,
                  Z = Z,
                  R = R,
                  index_obs = index_obs,
                  index_mis = index_mis,
                  N_obs = N_obs,
                  N_mis = N_mis)
  
  return(datlist)
}

