# Calculation of level of spatial variation in relationships
# By taxon, response variable, and predictor variable
# Approach: mean pairwise difference of coefficients among adjacent regions

# QDR NASAbioXgeo 19 Jun 2018

# Do with mean coefficients only, then later can get CIs.

# Load coefficients and adjacency matrix.
fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/modelfits' # Local
model_coef <- read.csv(file.path(fp, 'multivariate_spatial_coef.csv'), stringsAsFactors = FALSE)
load(file.path(fp,'fia_spatial_mm_dat_50k.RData'))

# Calculation of spatial variation of a coefficient

x <- subset(model_coef, model == 'full' & taxon == 'fia' & effect == 'coefficient' & response == 'alpha_richness' & parameter == 'elevation_5k_tri_50_mean' & stat == 'Estimate')

coefdiff <- as.matrix(dist(x$value)) # pairwise difference among coefficients
dimnames(coefdiff) <- list(x$region, x$region) # give region names
tnc_index <- match(x$region, dimnames(tnc_bin)[[1]]) # subset the TNC matrix by the 
tnc_used <- tnc_bin[tnc_index, tnc_index]

# Get only the entries that are neighbors and in the upper triangle of the matrix.
coefdiff[tnc_used == 0 | lower.tri(coefdiff)] <- NA
mean(na.omit(as.numeric(coefdiff)))

coef_diff <- function(x) {
  coefdiff <- as.matrix(dist(x$value)) # pairwise difference among coefficients
  dimnames(coefdiff) <- list(x$region, x$region) # give region names
  tnc_index <- match(x$region, dimnames(tnc_bin)[[1]]) # subset the TNC matrix by the 
  tnc_used <- tnc_bin[tnc_index, tnc_index]
  
  # Get only the entries that are neighbors and in the upper triangle of the matrix.
  coefdiff[tnc_used == 0 | lower.tri(coefdiff)] <- NA
  data.frame(coef_var = mean(na.omit(as.numeric(coefdiff))))
}

library(dplyr)

model_coef_var <- model_coef %>%
  filter(model == 'full', effect == 'coefficient', stat == 'Estimate') %>%
  group_by(taxon, response, parameter) %>%
  do(coef_diff(.)) %>%
  ungroup

write.csv(model_coef_var, file.path(fp, 'multivariate_spatial_coef_variation.csv'), row.names = FALSE)
