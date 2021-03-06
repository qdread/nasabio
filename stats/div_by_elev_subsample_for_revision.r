# Conceptual paper model fitting: FIA a,b,g diversity ~ elevation SD
# edited 10 Nov: add beta-diversity and also add gam fit with no subsample.
# edited 05 Dec: now elevation stats are correct, beta-diversity is based on pairwise and not multisite Soerensen.
# new version created 06 Dec: parallel so that many iterations can be run simultaneously on cluster
# Modified 28 Aug 2018: To comply with reviewer suggestions, change analysis to GLM with gamma distributions (and reduce sample number by a few)
# Modified 06 Sep 2018: Use standardized coefficients.
# Modified 28 Nov 2018: New input data (macroplot trees removed). This requires first getting rid of the non-PNW plots.

# Procedure in pseudocode:
# Load a,b,g diversity data and elevation diversity data
# Throw out all points within 100 km of Mexico and Canada.
# for iteration in 1:n_iter { 
#   Do iterative search to get ~20 points from the dataset (since that is as many 100 km radius ones as we can get)
#   for radius in 5, 10, 20, 50, 100 {
#     for y_variable in alpha, beta, gamma {
#       
#       Fit lm or betareg to those points, using the correct radius x-var and y-var
#       Save result to array
#       Get predicted values and save to array
# }}}

task <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
fpgeo <- '/mnt/research/nasabio/data/fia/geodiversity_CSVs'
fpbio <- '/mnt/research/nasabio/data/fia/biodiversity_CSVs/updated_nov2018'

# Load the elevational diversity and abg diversity data
ed <- read.csv(file.path(fpgeo, 'fia_usa_elev_only.csv'))
ad <- read.csv(file.path(fpbio, 'fiausa_natural_alpha.csv'))
bd <- read.csv(file.path(fpbio, 'fiausa_natural_betatd.csv'))
gd <- read.csv(file.path(fpbio, 'fiausa_natural_gamma.csv'))

# Load state codes
fia_statecodes <- read.csv('/mnt/research/nasabio/data/fia/fiastatecodes.csv')

radii <- c(5, 10, 20, 50, 100)
div_names <- c('alpha_diversity','beta_diversity','gamma_diversity')

mylib <- '/mnt/home/qdr/R/x86_64-pc-linux-gnu-library/3.5'
library(dplyr)
library(sp)
library(mgcv)
library(reghelper, lib.loc = mylib) # For standardized coefficients
library(betareg, lib.loc = mylib)

# Function for iterative search.
source('/mnt/research/nasabio/code/SRS_iterative.r')
# Functions for flagging edge plots
source('/mnt/research/nasabio/code/spatial_fns.r')

# Filter ed dataframe by state code to get only the PNW states, and only the non plantation plots.
pnw_states <- c(6, 41, 53) # Cal, Ore, and Wash.
ed <- ed %>%
  filter(PLT_CN %in% ad$PLT_CN) %>%
  left_join(fia_statecodes) %>%
  filter(STATECD %in% pnw_states) 

# Combine into a single data frame.
biogeo <- ed %>%
  dplyr::select(PLT_CN, STATECD, COUNTYCD, radius, sd) %>%
  filter(radius %in% radii) %>%
  rename(elevation_sd = sd) %>%
  left_join(ad %>% 
              dplyr::select(PLT_CN, radius, richness, shannon) %>% 
              rename(alpha_richness = richness, alpha_diversity = shannon) %>%
              filter(radius %in% radii)) %>%
  left_join(bd %>% 
              dplyr::select(PLT_CN, radius, beta_td_pairwise_pa, beta_td_pairwise) %>% 
              rename(beta_richness = beta_td_pairwise_pa, beta_diversity = beta_td_pairwise) %>%
              filter(radius %in% radii)) %>%
  left_join(gd %>% 
              dplyr::select(PLT_CN, radius, richness, shannon) %>% 
              rename(gamma_richness = richness, gamma_diversity = shannon) %>%
              filter(radius %in% radii))

# Add latitude and longitudes from unfuzzed (on local drive only)
fiacoords <- read.csv('/mnt/home/qdr/data/allfia.csv') %>% 
  rename(PLT_CN = CN) %>%
  left_join(fia_statecodes) %>%
  filter(STATECD %in% pnw_states)

# Convert to albers
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
fia_aea <- spTransform(SpatialPoints(coords=fiacoords[,c('ACTUAL_LON','ACTUAL_LAT')], proj4string = CRS('+proj=longlat')), CRSobj = CRS(aea_crs))

# Throw out the ones that are within 100 km of Mexico or Canada or the eastern edge
fia_coast100 <- flag_coast_plots(focal_points = fia_aea, radius = 100e3, focal_states = c('California','Oregon','Washington','Alaska'), border_states = c('Arizona','Idaho','Nevada'), border_countries = c('Canada', 'Mexico'))

fiacoords <- cbind(fiacoords, fia_coast100)

biogeo <- biogeo %>%
  left_join(fiacoords %>% dplyr::select(PLT_CN, ACTUAL_LAT, ACTUAL_LON, is_edge) %>% rename(lat = ACTUAL_LAT, lon = ACTUAL_LON)) %>%
  filter(!is_edge, complete.cases(.))
  
fia_aea_noedge <- fia_aea[!fiacoords$is_edge & fiacoords$PLT_CN %in% biogeo$PLT_CN]
fiacoords_noedge <- subset(fiacoords, !is_edge & PLT_CN %in% biogeo$PLT_CN)


# Calculate pairwise distances. (only takes ~1 min with brute force method)
fia_pnw_dist <- spDists(fia_aea_noedge, longlat = FALSE)

# Subsampling with glms ---------------------------------------------------

n_iter <- 1000 # Per job. Should take less than 4 hours to complete.
subsample_size <- 20
n_max <- 50
n_predict <- 50

# Get rid of the (very few) places where beta diversity is exactly 1
biogeo$beta_diversity[biogeo$beta_diversity == 1] <- NA

# Transform alpha and gamma diversity to "true" equivalents
biogeo <- mutate(biogeo,
                alpha_diversity = exp(alpha_diversity),
                gamma_diversity = exp(gamma_diversity))

# x-values for predicted y-values
xrange <- range(biogeo$elevation_sd, na.rm=TRUE)
newx <- round(seq(xrange[1], xrange[2], length.out = n_predict))

# R2 array with number of diversity types (a,b,g) by number of radii by number of iterations
r2_lm_array <- array(NA, dim = c(length(div_names), length(radii), n_iter))
# Added 15 Feb: Save slopes or coefficients of beta regression as well.
coef_array <- array(NA, dim = c(length(div_names), length(radii), n_iter))
# Predicted value array (n of diversity types by number of radii by number of iterations by number of predicted values)
pred_val_lm_array <- array(NA, dim = c(length(div_names), length(radii), n_iter, n_predict))

pb <- txtProgressBar(0, prod(dim(r2_lm_array)), style = 3)
ii <- 0
set.seed(task + 10101)

for (k in 1:n_iter) {
  
  # Subsample using the largest radius (100 km)
  sample_idx <- SRS_iterative(focal_points = fia_aea_noedge, dist_mat = fia_pnw_dist, radius = 100e3, n = n_max)
  
  # If too many are subsampled, take a random subset of the subsample.
  if (length(sample_idx) > subsample_size) {
    plots_in_sample <- fiacoords_noedge$PLT_CN[sample(sample_idx, subsample_size)]
  } else {
    plots_in_sample <- fiacoords_noedge$PLT_CN[sample_idx]
  }
  
  for (j in 1:length(radii)) {

    dat <- biogeo %>%
      filter(radius == radii[j], PLT_CN %in% plots_in_sample)
    
    for (i in 1:length(div_names)) {
      ii <- ii + 1
      setTxtProgressBar(pb, ii)
      
      # Fit lm if alpha or gamma, fit beta-regression (oddly enough) if beta
      if (div_names[i] == 'beta_diversity') {
        betareg_fit <- betareg(formula(paste(div_names[i], 'elevation_sd', sep = '~')), data = dat)
        pred_val_lm_array[i, j, k, ] <- predict(object = betareg_fit, newdata = data.frame(elevation_sd = newx))
        r2_lm_array[i, j, k] <- betareg_fit$pseudo.r.squared
		    coef_array[i, j, k] <- coef(betareg_fit)[2] * sd(dat$elevation_sd, na.rm = TRUE) # Manually standardized
      } else { 
        lm_fit <- glm(formula(paste(div_names[i], 'elevation_sd', sep = '~')), data = dat, family = Gamma(link = 'log'))
        pred_val_lm_array[i, j, k, ] <- predict.glm(object = lm_fit, newdata = data.frame(elevation_sd = newx))
        r2_lm_array[i, j, k] <- with(summary(lm_fit), 1 - deviance/null.deviance)
		    coef_array[i, j, k] <- coef(reghelper::beta(lm_fit))[2, 1] # Standardized using reghelper::beta()
      }

    }
  }
}

close(pb)

save(r2_lm_array, pred_val_lm_array, coef_array, file = paste0('/mnt/research/nasabio/temp/pnwfit/fit_', task, '.r'))
