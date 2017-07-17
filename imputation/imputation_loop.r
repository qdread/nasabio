# Create many missing datasets and run the rphylopars and mice imputations.
# Save the missing datasets so that we can pass them into the stan model.
# Calculate RMSEs for each trait separately, and for the entire model.

pct_missing <- 0.25 # Let's try this for now.
n_datasets <- 999 # If we cannot run 999 copies of the stan model, we will just use a subset of these.

# WORKFLOW
# 1. Generate missing datasets and save
# 2. Run the mice imputation n times and save the RMSE
# 3. Run the phylopars imputation n times and save the RMSE
# 4. Get RMSEs from STAN fits too.
# 5. Put RMSEs into a data frame and plot them for the "money" figure on the poster/in the paper.

# FUNCTION
# RMSE for each trait separately and overall

RMSE_eachtrait <- function(df) {
  RMSE_overall <- with(df, (sum(true_trait - imputed_trait)^2)/length(true_trait))
  require(dplyr)
  RMSE_bytrait <- df %>% group_by(trait_id) %>% summarize(RMSE = (sum(true_trait - imputed_trait)^2)/length(true_trait))
  result <- c(RMSE_overall, RMSE_bytrait$RMSE)
  names(result) <- c('overall', 'Bark.thickness', 'Wood.density', 'Specific.leaf.area', 'Plant.height', 'Plant.lifespan', 'Seed.dry.mass')
  result
}

# 1. generate missing datasets --------------------------------------------

source('~/GitHub/r_code/remove_values.R')
source('~/GitHub/nasabio/imputation/loadtraitsandpredictors.r')
source('~/GitHub/nasabio/imputation/makestandata.r')

# Log-transform trait data output by the above scripts (species_traits object with 64 species)

species_traits[,-1] <- log(species_traits[,-1])

# Generate missing values
set.seed(48824)

missing_datasets <- list()

for (i in 1:n_datasets) {
  missing_datasets[[i]] <- removeValues(species_traits, proportion = pct_missing)
}

# Save missing datasets on HPCC (commented out)
# save(missing_datasets, file = 'X:/data/fia/missing_datasets.R') # For future use.


# 2. Run mice imputation on missing datasets ------------------------------

mice_rmse <- function(dat, n_iter) {
  require(mice)
  init = mice(dat,maxit=50,method="pmm", print = FALSE)
  meth = init$method
  predM = init$predictorMatrix

  imputed <- mice(data = dat, method="pmm", m=n_iter, predictorMatrix = predM, print = FALSE) 
  
  # Extract values from the mice object and calculate standard errors for confidence intervals.
  imputed_raw <- imputed$imp[-1] # element 1 is removed because it was a column of species names.
  
  # Function to extract summary stats (row-wise mean and variance)
  get_mice_stats <- function(dat) {
    data.frame(mean = apply(dat, 1, mean),
               var = apply(dat, 1, var))
  }
  
  mice_summstats <- lapply(imputed_raw, get_mice_stats)
  
  
  # Put the imputed values into the missing data frame
  mice_complete <- complete(imputed)
  
  # Generate data frame with all true and imputed values side-by-side so that we can compare them.
  
  imputed_traits <- mice_complete[,-1]
  true_traits <- species_traits
  is_missing <- is.na(dat[,-c(1)])
  
  # The variances can be concatenated from the list of results
  mice_imputed_variances <- do.call(c, lapply(mice_summstats, '[', , 'var'))
  
  
  # Dataframe with comparison of values.
  comparisondf <- data.frame(imputed_trait =(imputed_traits[is_missing]),
                             imputed_variance = mice_imputed_variances,
                             true_trait = true_traits[, -1][is_missing],
                             trait_id = col(is_missing)[is_missing],  # Vector of column numbers with missing data 
                             species_id = row(is_missing)[is_missing])
  
  return(RMSE_eachtrait(comparisondf))
}



library(pbapply)
mice_RMSEs <- pblapply(missing_datasets, mice_rmse, n_iter = 100) # Each of the 999 datasets will be imputed 100 times. Run time ~ 2 hours

save(mice_RMSEs, file = 'insert path')

# 3. Run Rphylopars imputation on missing datasets ------------------------

rphylopars_rmse <- function(dat, phy) {
  require(Rphylopars)
  require(ape)
  
  # make tree multi2di (collapse or resolve multichotomies in phylogenetic trees)
  phy_corrected <- multi2di(phy)
  
  phy_OU <- phylopars(dat, phy_corrected, model='OU')
  
  # Generate the comparison dataframe.
  n_spp <- length(unique(dat$species))
  
  imputed_traits <- phy_OU$anc_recon[1:n_spp,] # Get only rows from existing species, not the ancestral nodes species
  imputed_variances <- phy_OU$anc_var[order(rownames(phy_OU$anc_var[1:n_spp,])),] # Order the species names alphabetically
  ordered_trait_data <- phy_OU$trait_data[order(phy_OU$trait_data$species),] # Order trait data species names alphabetically
  is_missing <- is.na(ordered_trait_data[,-c(1)])
  true_traits <- species_traits
  
  comparisondf <- data.frame(imputed_trait =(imputed_traits[is_missing]),
                             imputed_variance = (imputed_variances[is_missing]),
                             true_trait = true_traits[, -1][is_missing],
                             trait_id = col(is_missing)[is_missing],  # Vector of column numbers with missing data 
                             species_id = row(is_missing)[is_missing])
  
  return(RMSE_eachtrait(comparisondf))
  
}

rphylo_RMSEs <- pblapply(missing_datasets, rphylopars_rmse, phy = species_phylogeny)

save(rphylo_RMSEs, file = 'insert path')

# 4. load stan fits and get rmses from them -------------------------------

# Make stan rdumps for each of the missing datasets and write them to .R files on the hpcc
library(rstan)

for (i in 1:n_datasets) {
  missing_list_i <- make_standatalist_missing(traitdat = missing_datasets[[i]], predictors = species_covariates, phy = species_phylogeny, evolution_model = 'ou')
  with(missing_list_i, stan_rdump(names(missing_list_i), file = paste0('X:/data/fia/staninputs/trait_data_list_', i, '.R')))
}



# 5. Combine RMSEs into data frame and plot -------------------------------

mice_RMSEs <- do.call(rbind, mice_RMSEs)
rphylo_RMSEs <- do.call(rbind, rphylo_RMSEs)
# stan rmses will be added here when ready.

library(reshape2)
mice_RMSEs_long <- melt(mice_RMSEs, value.name = 'RMSE', variable.name = 'trait')
rphylo_RMSEs_long <- melt(rphylo_RMSEs, value.name = 'RMSE', variable.name = 'trait')
RMSE_data <- rbind(data.frame(method = 'MICE', mice_RMSEs_long),
                   data.frame(method = 'Rphylopars', rphylo_RMSEs_long))

write.csv(RMSE_data, file = '...csv', row.names = FALSE)

library(cowplot)

# For now, plot just the overall results. Might also be worth looking at the plots for the individual trait results as well.
ggplot(subset(RMSE_data, trait == 'overall'), aes(x = method, y = RMSE)) +
  geom_boxplot() +
  labs(x = 'Imputation method', y = 'Root mean squared error')