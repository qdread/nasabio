# MICE imputation with code to calculate mean and variance of imputed values
# and format them to pass into ggplot2
# QDR modified JJ's code on 17 July 2017

library(mice)
training <- read.csv('X:/jay/tree_data/stevenstraits/trait_stevens_training.csv',head = TRUE, sep = ',', check.names = FALSE, stringsAsFactors = FALSE)

# Take log
miceLogTraining <- training
miceLogTraining[,-1] <- log(miceLogTraining[,-1])

# Remove values at random
source('~/GitHub/r_code/remove_values.R')
set.seed(1719)
miceLogTraining <- removeValues(miceLogTraining, proportion = 0.2)

init = mice(miceLogTraining,maxit=50,method="pmm")
meth = init$method
predM = init$predictorMatrix

# Use 999 iterations to get robust distribution around mean of missing value
# In this test version, just do 10 iterations.
imputed <- mice(data = miceLogTraining, method="pmm", m=10, predictorMatrix = predM) 

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


#imputed_traits <- imputed # Get only rows from existing species, not the ancestral nodes species
#imputed_variances <- phy_training_OU$anc_var[order(rownames(phy_training_OU$anc_var[1:83,])),] # Order the species names alphabetically
#ordered_trait_data <- phy_training_OU$trait_data[order(phy_training_OU$trait_data$species),] # Order trait data species names alphabetically
imputed_traits <- mice_complete[,-1]
true_traits <- training
is_missing <- is.na(miceLogTraining[,-c(1)])

# The variances can be concatenated from the list of results
mice_imputed_variances <- do.call(c, lapply(mice_summstats, '[', , 'var'))


#plot(x = training[,-1], y = imputed[,-1])

# Dataframe with comparison of values.
comparisondf <- data.frame(imputed_trait =(imputed_traits[is_missing]),
                           imputed_variance = mice_imputed_variances,
                           true_trait = log(true_traits[-1][is_missing]),
                           trait_id = col(is_missing)[is_missing],  # Vector of column numbers with missing data 
                           species_id = row(is_missing)[is_missing])

traitNames <- c('1'='Bark thickness', '2'='Wood density', '3'='SLA','4' = 'Plant height','5' = 'Plant lifespan','6' = 'Seed dry mass')

comparisondf$speciesName <- training$Scientific_Name[comparisondf$species_id]

# Plot results. The red point is the "true" value. It should fall inside the confidence interval around the black point. It seems to perform pretty well for the fake data.
# Confidence intervals are based on assumption of normal distribution.
library(ggplot2)
library(cowplot)

all_traits <- ggplot(comparisondf, aes(x = speciesName)) + facet_wrap(~ trait_id, scales = 'free',labeller = labeller(trait_id = traitNames)) +
  geom_pointrange(aes(y = imputed_trait, 
                      ymin = imputed_trait - 1.96 * sqrt(imputed_variance), 
                      ymax = imputed_trait + 1.96 * sqrt(imputed_variance))) +
  geom_point(aes(y = true_trait), color = 'red', shape = 19) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5))+
  ggtitle("Imputations with 95% CI at 20% Missing Values using Training Dataset\n") +
  labs( x = "Species names", y = "Trait values") 
