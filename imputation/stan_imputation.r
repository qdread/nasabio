#####################################################
# Step 1: compile stan model

library(rstan)

trait_model <- stanc(file = 'imputation/phylo_spatial_trait.stan') # Read from GitHub.

#####################################################
# Step 2: Create list containing m,n,N,p,Y,X,Z,and R for model fitting

# environmental predictors, phylogenetic VCV matrix, trait data, and spatial locations

# R: phylogenetic VCV matrix. Can be either Brownian or Ornstein.
library(ape)
allfia_phylogeny <- read.tree('X:/data/fia/tree_all_final_031716.txt')

vcv_brownian <- vcv(allfia_phylogeny)

# Function for vcv under OU model is in PIGShift
library(PIGShift)
vcv_ou <- OU.vcv(allfia_phylogeny, theta = 1)

# For now, create environmental predictors and traits from species means (mean locations)

try_indiv_traits <- read.csv('X:/data/fia/fia_try_04jul/try_trait_byobs_all.csv', stringsAsFactors = FALSE)


# Very fake data
fakespmeans <- cbind(tr1 = c(5,4,6), tr2 = c(1,1,5), tr3 = c(4,3,2))
npersp <- 10
traitdat <- data.frame(species = rep(c('sp1','sp2','sp3'), each=npersp*3), 
                      trait = rep(c('tr1','tr2','tr3'), each=npersp),
                      value = c(rnorm(30, rep(fakespmeans[,'tr1'], each=10), 1),
                                rnorm(30, rep(fakespmeans[,'tr2'], each=10), 1),
                                rnorm(30, rep(fakespmeans[,'tr3'], each=10), 1)))

# Z is design matrix with dimensions (n trait * n indiv) x (n trait * n species) to make sure each individual maps to the right species
Z <- matrix(0, nrow = nrow(traitdat), ncol = length(unique(traitdat$trait)) * length(unique(traitdat$species)))

# For each row in traitdat, put a 1 in the correct row and column of Z.
species_trait_combo <- with(traitdat, paste(species,trait))
col_idx <- match(species_trait_combo, unique(species_trait_combo))

for (i in 1:length(col_idx)) Z[i, col_idx[i]] <- 1

trait_data_list <- list(m = length(unique(traitdat$trait)),
                        n = length(unique(traitdat$species)),
                        N = nrow(traitdat),
                        p = ncol(X) + 1,
                        X = as.matrix(X),
                        Y = as.matrix(Y),
                        Z = Z,
                        R = as.matrix(R))

######################################################
# Step 3: Set model fitting options

######################################################
# Step 4: Fit model

######################################################
# Step 5: Examine diagnostic plots and metrics to see if the model converged