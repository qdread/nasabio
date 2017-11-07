# Conceptual paper model fitting: FIA a,b,g diversity ~ elevation SD

# Procedure in pseudocode:
# Load a,b,g diversity data and elevation diversity data
# Throw out all points within 100 km of Mexico and Canada.
# for iteration in 1:999 { 
#   Do iterative search to get ~20 points from the dataset (since that is as many 100 km radius ones as we can get)
#   for radius in 5, 10, 20, 50, 100 {
#     for y_variable in alpha, beta, gamma {
#       
#       Fit GAM to those points, using the correct radius x-var and y-var
#       Save result to some sort of array
# }}}


# NOTE: All needed CSVs are on the hpcc but I have downloaded them locally because it is faster. 
# Change file path to the second file path to get files from hpcc.

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/fia_unfuzzed'
#fp <- '/mnt/research/nasabio/data/fia'

# Load the elevational diversity and abg diversity data
ed <- read.csv(file.path(fp, 'fia_elev_stats_unfuzzed.csv'))
ad <- read.csv(file.path(fp, 'fia_alpha.csv'))
gd <- read.csv(file.path(fp, 'fia_gammadiv.csv'))

radii <- c(5, 10, 20, 50, 100)
div_names <- c('alpha_richness','beta_richness','gamma_richness')

library(dplyr)
library(sp)
library(mgcv)

# Function for iterative search.
source('~/GitHub/nasabio/methods/SRS_iterative.r')
# Functions for flagging edge plots
source('~/GitHub/nasabio/methods/spatial_fns.r')

# Combine into a single data frame.
biogeo <- ed %>%
  dplyr::select(PLT_CN, radius, sd) %>%
  filter(radius %in% radii) %>%
  rename(elevation_sd = sd) %>%
  left_join(ad %>% dplyr::select(PLT_CN, radius, richness) %>% rename(alpha_richness = richness)) %>%
  left_join(gd %>% dplyr::select(PLT_CN, radius, richness) %>% rename(gamma_richness = richness))

# Add latitude and longitudes from unfuzzed (on local drive only)
fiacoords <- read.csv('~/FIA/pnw.csv') %>% filter(complete.cases(.))

# Convert to albers
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
fia_aea <- spTransform(SpatialPoints(coords=fiacoords[,c('ACTUAL_LON','ACTUAL_LAT')], proj4string = CRS('+proj=longlat')), CRSobj = CRS(aea_crs))

# Throw out the ones that are within 100 km of Mexico or Canada or the eastern edge
fia_coast100 <- flag_coast_plots(focal_points = fia_aea, radius = 100e3, focal_states = c('California','Oregon','Washington','Alaska'), border_states = c('Arizona','Idaho','Nevada'), border_countries = c('Canada', 'Mexico'))

fiacoords <- cbind(fiacoords, fia_coast100)

biogeo <- biogeo %>%
  left_join(fiacoords %>% dplyr::select(CN, ACTUAL_LAT, ACTUAL_LON, is_edge) %>% rename(PLT_CN = CN, lat = ACTUAL_LAT, lon = ACTUAL_LON)) %>%
  filter(!is_edge, complete.cases(.))
  
fia_aea_noedge <- fia_aea[!fiacoords$is_edge & fiacoords$CN %in% biogeo$PLT_CN]
fiacoords_noedge <- subset(fiacoords, !is_edge & CN %in% biogeo$PLT_CN)


# Calculate pairwise distances. (only takes ~1 min with brute force method)
fia_pnw_dist <- spDists(fia_aea_noedge, longlat = FALSE)


# Subsampling with gams ---------------------------------------------------

n_iter <- 999
subsample_size <- 20
n_max <- 50

# R2 of gam. array with number of diversity types (a,b,g) by number of radii by number of iterations
r2_array <- array(NA, dim = c(length(div_names), length(radii), n_iter))
r2_lm_array <- array(NA, dim = c(length(div_names), length(radii), n_iter))

# x-values for predicted y-values
xrange <- range(biogeo$elevation_sd, na.rm=TRUE)
newx <- round(seq(xrange[1], xrange[2], length.out = 20))

# Predicted value array
pred_val_array <- array(NA, dim = c(length(div_names), length(radii), n_iter, length(newx)))
pred_val_lm_array <- array(NA, dim = c(length(div_names), length(radii), n_iter, length(newx)))

pb <- txtProgressBar(0, prod(dim(r2_array)), style = 3)
ii <- 0
set.seed(517)

for (k in 1:n_iter) {
  
  # Subsample using the largest radius (100 km)
  sample_idx <- SRS_iterative(focal_points = fia_aea_noedge, dist_mat = fia_pnw_dist, radius = 100e3, n = n_max)
  
  # If too many are subsampled, take a random subset of the subsample.
  if (length(sample_idx) > subsample_size) {
    plots_in_sample <- fiacoords_noedge$CN[sample(sample_idx, subsample_size)]
  } else {
    plots_in_sample <- fiacoords_noedge$CN[sample_idx]
  }
  
  for (j in 1:length(radii)) {

    dat <- biogeo %>%
      filter(radius == radii[j], PLT_CN %in% plots_in_sample)
    
    for (i in 1:length(div_names)) {
      ii <- ii + 1
      setTxtProgressBar(pb, ii)
      # Fit GAM.
      gam_fit <- gam(formula(paste(div_names[i], 'elevation_sd', sep = '~')), data = dat)
      lm_fit <- lm(formula(paste(div_names[i], 'elevation_sd', sep = '~')), data = dat)
      
      # If possible, get predicted y-values for the gam. This can be used to make an interval on the figures.
      pred_val_array[i, j, k, ] <- predict.gam(object = gam_fit, newdata = data.frame(elevation_sd = newx))
      pred_val_lm_array[i, j, k, ] <- predict.lm(object = lm_fit, newdata = data.frame(elevation_sd = newx))
      
      # Save r-squared.
      r2_array[i, j, k] <- summary(gam_fit)$r.sq
      r2_lm_array[i, j, k] <- summary(lm_fit)$r.sq
    }
  }
}

close(pb)

# Convert array to data frame
library(reshape2)
dimnames(r2_array) <- list(div_names, radii, NULL)
dimnames(r2_lm_array) <- list(div_names, radii, NULL)
dimnames(pred_val_array) <- list(div_names, radii, NULL, NULL)
dimnames(pred_val_lm_array) <- list(div_names, radii, NULL, NULL)


r2_df <- melt(r2_array, varnames = c('diversity_type', 'radius', 'iteration'))
pred_val_df <- melt(pred_val_array, varnames = c('diversity_type', 'radius', 'iteration', 'x'))
pred_val_df$x <- newx[pred_val_df$x]

r2_lm_df <- melt(r2_lm_array, varnames = c('diversity_type', 'radius', 'iteration'))
pred_val_lm_df <- melt(pred_val_lm_array, varnames = c('diversity_type', 'radius', 'iteration', 'x'))
pred_val_lm_df$x <- newx[pred_val_lm_df$x]

# Plot
library(cowplot)

### GAMs
# Variation explained
ggplot(r2_df, aes(x = radius, y = value, group = interaction(radius, diversity_type), fill = diversity_type)) +
  geom_boxplot()

# Predicted value
ggplot(pred_val_df, aes(x = x, y = value)) +
  facet_grid(diversity_type ~ radius, scales = 'free') +
  geom_line(aes(group = iteration))

ggplot(pred_val_df, aes(x = x, y = value)) +
  facet_grid(diversity_type ~ radius, scales = 'free') +
  stat_summary(geom = 'pointrange')

# Try to create a confidence ribbon around the predictive interval so that it looks smoother. Can also plot the points on that.
pred_val_quant <- pred_val_df %>%
  group_by(diversity_type, radius, x) %>%
  summarize(pred_y = quantile(value, probs = 0.5),
            pred_y_min = quantile(value, probs = 0.025),
            pred_y_max = quantile(value, probs = 0.975)) %>%
  arrange(diversity_type, radius, x)

ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free') +
  geom_point(aes(x = elevation_sd, y = alpha_richness), alpha = 0.05) +
  geom_ribbon(data = pred_val_quant %>% filter(diversity_type == 'alpha_richness'), 
              aes(x = x, ymin = pred_y_min, ymax = pred_y_max), fill = 'blue', alpha = 0.3) +
  geom_line(data = pred_val_quant %>% filter(diversity_type == 'alpha_richness'), 
            aes(x = x, y = pred_y), color = 'blue', size = 1.5)

ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free') +
  geom_point(aes(x = elevation_sd, y = gamma_richness), alpha = 0.05) +
  geom_ribbon(data = pred_val_quant %>% filter(diversity_type == 'gamma_richness'), 
              aes(x = x, ymin = pred_y_min, ymax = pred_y_max), fill = 'blue', alpha = 0.3) +
  geom_line(data = pred_val_quant %>% filter(diversity_type == 'gamma_richness'), 
            aes(x = x, y = pred_y), color = 'blue', size = 1.5)

### LMs

# Variation explained
ggplot(r2_lm_df, aes(x = radius, y = value, group = interaction(radius, diversity_type), fill = diversity_type)) +
  geom_boxplot()

# Predicted values
pred_val_lm_quant <- pred_val_lm_df %>%
  group_by(diversity_type, radius, x) %>%
  summarize(pred_y = quantile(value, probs = 0.5),
            pred_y_min = quantile(value, probs = 0.025),
            pred_y_max = quantile(value, probs = 0.975)) %>%
  arrange(diversity_type, radius, x)

ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free') +
  geom_point(aes(x = elevation_sd, y = alpha_richness), alpha = 0.05) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'alpha_richness'), 
              aes(x = x, ymin = pred_y_min, ymax = pred_y_max), fill = 'blue', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'alpha_richness'), 
            aes(x = x, y = pred_y), color = 'blue', size = 1.5)

ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free') +
  geom_point(aes(x = elevation_sd, y = gamma_richness), alpha = 0.05) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'gamma_richness'), 
              aes(x = x, ymin = pred_y_min, ymax = pred_y_max), fill = 'blue', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'gamma_richness'), 
            aes(x = x, y = pred_y), color = 'blue', size = 1.5)
