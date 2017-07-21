# Test stan summary output

fp <- '/mnt/research/nasabio/data/fia/stanoutputs/shorter'

RMSE_eachtrait <- function(df) {
  RMSE_overall <- with(df, (sum(true_trait - imputed_trait)^2)/length(true_trait))
  require(dplyr)
  RMSE_bytrait <- df %>% group_by(trait_id) %>% summarize(RMSE = (sum(true_trait - imputed_trait)^2)/length(true_trait))
  result <- c(RMSE_overall, RMSE_bytrait$RMSE)
  names(result) <- c('overall', 'Bark.thickness', 'Wood.density', 'Specific.leaf.area', 'Plant.height', 'Plant.lifespan', 'Seed.dry.mass')
  result
}

# environmental predictors, phylogenetic VCV matrix, trait data, and spatial locations
library(dplyr)

# For now, create environmental predictors and traits from species means (mean locations)

# Use Stevens traits, and get median locations of each species from FIA dataset.
# Use the environmental predictors from FIA.

trait_training <- read.csv('/mnt/research/nasabio/jay/tree_data/stevenstraits/trait_stevens_training.csv', stringsAsFactors = FALSE)

# Get median locations for each species from TRY
# All covariates for those locations have already been extracted from different GIS layers.
try_covariates <- read.csv('/mnt/research/nasabio/data/fia/imputation_covariates.csv', stringsAsFactors = FALSE)

# Get only a subset of the covariates: elevation, minimum temperature, temperature seasonality, minimum precipitation, precipitation seasonality, lai, npp
try_covariates_reduced <- try_covariates %>%
  select(AccSpeciesName, lat_aea, lon_aea, elevation, bio6, bio4, bio14, bio15, lai, npp) %>%
  rename(temp_min_monthly = bio6, temp_seasonality = bio4, precip_min_monthly = bio14, precip_seasonality = bio15) %>%
  filter(complete.cases(.))


library(ape)
allfia_phylogeny <- read.tree('/mnt/research/nasabio/data/fia/tree_all_final_031716.txt')

# Pare down the trait list and phylogeny to species that we have valid locations for. (64 spp)

species_list <- Reduce('intersect', list(try_covariates_reduced$AccSpeciesName, allfia_phylogeny$tip.label, trait_training$Scientific_Name))

species_covariates <- subset(try_covariates_reduced, AccSpeciesName %in% species_list)
species_phylogeny <- drop.tip(allfia_phylogeny, tip = allfia_phylogeny$tip.label[!allfia_phylogeny$tip.label %in% species_list])
species_traits <- subset(trait_training, Scientific_Name %in% species_list)

# Sort them all so that they are in the same order
species_covariates <- species_covariates %>%
  rename(species = AccSpeciesName) %>%
  arrange(species)

species_traits <- species_traits %>%
  rename(species = Scientific_Name) %>%
  arrange(species)

# Mean-center and standardize the predictor variables.
species_covariates[,-1] <- sapply(species_covariates[,-1], scale)

# Get rid of everything except temperature, precipitation, and productivity
species_covariates <- species_covariates[,c('species','temp_min_monthly','precip_min_monthly','npp')]

species_traits[,-1] <- log(species_traits[,-1])

load(file = '/mnt/research/nasabio/data/fia/missing_datasets.R')

# Loop through the list of loaded stan summaries and output a comparisondf for each one.

stan_comparisondfs <- list()
stan_RMSEs <- list()

summary_files <- dir(fp, pattern = 'summary*')
summary_numbers <- as.numeric(stringr::str_extract(summary_files, '[0-9]+'))  

for (i in summary_files) {

summary_number <- as.numeric(stringr::str_extract(i, '[0-9]+'))  
summary_i <- read.table(file.path(fp, i), skip=6, nrow=2000)
  
true_traits <- species_traits
is_missing <- is.na(missing_datasets[[summary_number]][,-c(1)])

y_miss <- summary_i[grepl('Y_mis', row.names(summary_i)), ]

imputed_traits <- y_miss$X50.
imputed_ci_min <- y_miss$X5.
imputed_ci_max <- y_miss$X95.

# Dataframe with comparison of values.

comparisondf <- data.frame(imputed_trait =imputed_traits,
                           imputed_ci_min = imputed_ci_min,
                           imputed_ci_max = imputed_ci_max,
                           true_trait = true_traits[, -1][is_missing],
                           trait_id = col(is_missing)[is_missing],  # Vector of column numbers with missing data 
                           species_id = row(is_missing)[is_missing])

stan_comparisondfs[[length(stan_comparisondfs) + 1]] <- comparisondf
stan_RMSEs[[length(stan_RMSEs) + 1]] <- RMSE_eachtrait(comparisondf)

}

# Combine the stan rmses together

stan_RMSEs <- as.data.frame(do.call(rbind, stan_RMSEs))

save(stan_RMSEs, summary_numbers, file = '/mnt/research/nasabio/jay/stan_rmse_object.r')

load('/mnt/research/nasabio/jay/stan_rmse_object.r')

mice_RMSEs <- mice_RMSEs[summary_numbers, ]
rphylo_RMSEs <- rphylo_RMSEs[summary_numbers, ]

library(reshape2)
mice_RMSEs_long <- melt(mice_RMSEs, value.name = 'RMSE', variable.name = 'trait')
rphylo_RMSEs_long <- melt(rphylo_RMSEs, value.name = 'RMSE', variable.name = 'trait')
stan_RMSEs_long <- melt(stan_RMSEs, value.name = 'RMSE', variable.name = 'trait')

RMSE_data <- rbind(data.frame(method = 'MICE', mice_RMSEs_long),
                   data.frame(method = 'Rphylopars', rphylo_RMSEs_long),
                   data.frame(method = 'Hierarchical', stan_RMSEs_long))

