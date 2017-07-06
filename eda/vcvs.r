# Variance-covariance matrix for FIA phylogeny assuming different models of evolution: Brownian and O-U

# Phylogeny
library(ape)
allfia_phylogeny <- read.tree('C:/Users/Q/Dropbox/projects/nasabiodiv/allfiaphylogeny/tree_all_final_031716.txt')

vcv_brownian <- vcv(allfia_phylogeny)

# Function for vcv under OU model is in PIGShift
library(PIGShift)
vcv_ou <- OU.vcv(allfia_phylogeny, theta = 1)

# Compute eigenvalues
any(eigen(vcv_brownian)$values < 0)
any(eigen(vcv_ou)$values < 0)
# Matrices are positive definite.