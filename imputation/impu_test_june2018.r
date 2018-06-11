# Simple imputation test.
# Compare MICE built in method with own code written in Stan.

traits <- read.csv('~/R/trait_stevens_training.csv', stringsAsFactors = FALSE)
traits_cat <- read.csv('~/R/trait_stevens_categorical.csv', stringsAsFactors = FALSE)
traits_con <- read.csv('~/R/trait_stevens_continuous.csv', stringsAsFactors = FALSE)

dimnames(traits)[[1]] <- traits$Scientific_Name
traits <- traits[,-1]

set.seed(313)
pct_miss <- 0.2

n_miss <- round(pct_miss * prod(dim(traits)))

traits_MCAR <- as.matrix(traits)
traits_MCAR[sample(prod(dim(traits_MCAR)), size = n_miss)] <- NA
traits_MCAR <- as.data.frame(traits_MCAR)

library(mice)

traits_MICE_default <- mice(traits_MCAR, m = 5)
traits_MICE_bayesnorm <- mice(traits_MCAR, m = 5, method = 'norm')
