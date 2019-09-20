# Exploratory script to see whether there are any issues with too many correlated traits in the FD index.

birdtrait <- read.csv('~/Dropbox/projects/nasabiodiv/birdtraitmerged.csv')

birdtrait_imputed <- read.csv('~/Dropbox/projects/nasabiodiv/bbs_traits_imputed.csv')

library(FD)

birdtrait[birdtrait == -999] <- NA
nocturnalbirdAOUs <- birdtrait$AOU[birdtrait$Nocturnal == 1]
birdtrait_diurnal <- subset(birdtrait, Nocturnal != 1 | is.na(Nocturnal))
traitnames <- names(birdtrait)[c(15:24, 29:36, 46:50, 53, 55:59)]


# Explore paired correlations
tr <- birdtrait_diurnal[, traitnames]
tr_imp <- birdtrait_imputed[, traitnames]

tr_cor <- cor(tr_imp, use = 'pairwise.complete.obs')

# Generate heat map from fn in corheatmap.r
cormat_heatmap(round(tr_cor, 2))

# Identify the traits most correlated with adult body mass
cors_sorted <- sort(abs(tr_cor['adult_body_mass_g', ]))

# Compare the distance matrix calculated from ordination from the distance matrix calculated from the raw traits
# First get rid of the rows with a couple NAs
tr_imp_comp <- tr_imp[complete.cases(tr_imp), ]
tr_pca <- prcomp(tr_imp_comp, center = TRUE, scale = TRUE)

# Calculate distance matrices
dist_full <- gowdis(tr_imp_comp)
dist_pca5 <- gowdis(tr_pca$x[, 1:5])
dist_pca10 <- gowdis(tr_pca$x[, 1:10])

mantel(dist_full, dist_pca5, permutations = 999)
mantel(dist_full, dist_pca10, permutations = 999)

# Also try one that gets rid of the worst offending correlated variables
#usdm::vifcor(tr_imp_comp)
worst <- c('egg_mass_g', 'male_maturity_d', 'longevity_y', 'birth_or_hatching_weight_g')
tr_imp_reduced <- tr_imp_comp[, !names(tr_imp_comp) %in% worst]
cormat_heatmap(round(cor(tr_imp_reduced), 2))

dist_reduced <- gowdis(tr_imp_reduced)

mantel(dist_full, dist_reduced, permutations = 999)
