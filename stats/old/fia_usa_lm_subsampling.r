# Conceptual paper model fitting: FIA a,b,g diversity ~ elevation SD
# edited 10 Nov: add beta-diversity and also add gam fit with no subsample.
# edited 05 Dec: now elevation stats are correct, beta-diversity is based on pairwise and not multisite Soerensen.
# forked 09 Jan: same procedure, but use entire USA data, not just Pacific Northwest

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
ed <- read.csv(file.path(fp, 'fia_usa_elev_only.csv'))
ad <- read.csv(file.path(fp, 'fiausa_alpha.csv'))
bd <- read.csv(file.path(fp, 'fiausa_betatd.csv'))
gd <- read.csv(file.path(fp, 'fiausa_gamma.csv'))

# Correct errors in ed
ed <- filter(ed, min > -200, max < 5000)

radii <- c(5, 10, 20, 50, 100)
div_names <- c('alpha_diversity','beta_diversity','gamma_diversity')

library(dplyr)
library(sp)
library(mgcv)

# Function for iterative search.
source('~/GitHub/nasabio/stats/SRS_iterative.r')
# Functions for flagging edge plots
source('~/GitHub/nasabio/stats/spatial_fns.r')

# Combine into a single data frame.
biogeo <- ed %>%
  dplyr::select(PLT_CN, radius, sd) %>%
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
fiacoords <- read.csv('~/FIA/FIA10nov/allfia.csv') %>% filter(complete.cases(.))

# Convert to albers
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
fia_aea <- spTransform(SpatialPoints(coords=fiacoords[,c('ACTUAL_LON','ACTUAL_LAT')], proj4string = CRS('+proj=longlat')), CRSobj = CRS(aea_crs))

# Throw out the ones that are within 100 km of Mexico or Canada or the eastern edge
fia_coast100 <- flag_coast_plots(focal_points = fia_aea, radius = 100e3, border_countries = c('Canada', 'Mexico'))

fiacoords <- cbind(fiacoords, fia_coast100)

biogeo <- biogeo %>%
  left_join(fiacoords %>% dplyr::select(CN, ACTUAL_LAT, ACTUAL_LON, is_edge) %>% rename(PLT_CN = CN, lat = ACTUAL_LAT, lon = ACTUAL_LON)) %>%
  filter(!is_edge, complete.cases(.))
  
fia_aea_noedge <- fia_aea[!fiacoords$is_edge & fiacoords$CN %in% biogeo$PLT_CN]
fiacoords_noedge <- subset(fiacoords, !is_edge & CN %in% biogeo$PLT_CN)

# Get rid of the (very few) places where beta diversity is exactly 1
biogeo$beta_diversity[biogeo$beta_diversity == 1] <- NA

# Transform alpha and gamma diversity to "true" equivalents
biogeo <- mutate(biogeo,
                 alpha_diversity = exp(alpha_diversity),
                 gamma_diversity = exp(gamma_diversity))


# Fit LMs with no subsampling --------------------------------------------

# This is just to return the R-squared for the LMs. Added 10 Nov.

biogeo %>% 
  group_by(radius) %>%
  summarize(r2_alpharich = summary(lm(alpha_richness ~ elevation_sd))$r.sq,
            r2_betarich = summary(lm(beta_richness ~ elevation_sd))$r.sq,
            r2_gammarich = summary(lm(gamma_richness ~ elevation_sd))$r.sq,
            r2_alphadiv = summary(lm(alpha_diversity ~ elevation_sd))$r.sq,
            r2_betadiv = summary(lm(beta_richness ~ elevation_sd))$r.sq,
            r2_gammadiv = summary(lm(gamma_diversity ~ elevation_sd))$r.sq)

# Subsampling with LMs ---------------------------------------------------

n_iter <- 999 # Increase later.
subsample_size <- 150
n_max <- 200

library(betareg)


# x-values for predicted y-values
xrange <- range(biogeo$elevation_sd, na.rm=TRUE)
newx <- round(seq(xrange[1], xrange[2], length.out = 50))

# R2 array with number of diversity types (a,b,g) by number of radii by number of iterations
#r2_array <- array(NA, dim = c(length(div_names), length(radii), n_iter))
# Predicted value array
#pred_val_array <- array(NA, dim = c(length(div_names), length(radii), n_iter, length(newx)))

r2_lm_array <- array(NA, dim = c(length(div_names), length(radii), n_iter))
pred_val_lm_array <- array(NA, dim = c(length(div_names), length(radii), n_iter, length(newx)))

pb <- txtProgressBar(0, prod(dim(r2_lm_array)), style = 3)
ii <- 0
set.seed(517)

for (k in 1:n_iter) {
  
  # Subsample using the largest radius (100 km)
  sample_idx <- SRS_iterative_N1(focal_points = fia_aea_noedge, radius = 100e3, n = n_max)
  
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
      #gam_fit <- gam(formula(paste(div_names[i], 'elevation_sd', sep = '~')), data = dat)
      # If possible, get predicted y-values for the gam. This can be used to make an interval on the figures.
      #pred_val_array[i, j, k, ] <- predict.gam(object = gam_fit, newdata = data.frame(elevation_sd = newx))
      
      # Fit lm if alpha or gamma, fit beta-regression (oddly enough) if beta
      if (div_names[i] == 'beta_diversity') {
        betareg_fit <- betareg(formula(paste(div_names[i], 'elevation_sd', sep = '~')), data = dat)
        pred_val_lm_array[i, j, k, ] <- predict(object = betareg_fit, newdata = data.frame(elevation_sd = newx))
        r2_lm_array[i, j, k] <- betareg_fit$pseudo.r.squared
      } else { 
        lm_fit <- lm(formula(paste(div_names[i], 'elevation_sd', sep = '~')), data = dat)
        pred_val_lm_array[i, j, k, ] <- predict.lm(object = lm_fit, newdata = data.frame(elevation_sd = newx))
        r2_lm_array[i, j, k] <- summary(lm_fit)$r.sq
      }

    }
  }
}

close(pb)

###############################################
# 1e6 iterations were run on the cluster. Load them there, calculate summary information, and save to make plots locally.
# Edit 16 Feb 2018: Just use 1e5 iterations, with slopes saved.

radii <- c(5, 10, 20, 50, 100)
div_names <- c('alpha_diversity','beta_diversity','gamma_diversity')

# x-values for predicted y-values
xrange <- c(2, 1211)
newx <- round(seq(xrange[1], xrange[2], length.out = 50))

fp <- '/mnt/research/nasabio/data/fia/modelfits'

library(dplyr)

pred_list <- list()
r2_list <- list()
slope_list <- list()

fitnames <- dir(fp, pattern='fit_')

for (i in 1:10) {
  load(file.path(fp, fitnames[i]))
  pred_list[[i]] <- pred_val_lm_array
  r2_list[[i]] <- r2_lm_array
  slope_list[[i]] <- coef_array
  print(i)
}

library(abind)
pred_val_lm_array <- do.call('abind', c(pred_list, along = 3)) # very big.
r2_lm_array <- do.call('abind', c(r2_list, along = 3))
coef_array <- do.call('abind', c(slope_list, along=3))

# Convert array to data frame
library(reshape2)

#dimnames(r2_array) <- list(div_names, radii, NULL)
#dimnames(pred_val_array) <- list(div_names, radii, NULL, NULL)
#r2_df <- melt(r2_array, varnames = c('diversity_type', 'radius', 'iteration'))
#pred_val_df <- melt(pred_val_array, varnames = c('diversity_type', 'radius', 'iteration', 'x'))
#pred_val_df$x <- newx[pred_val_df$x]

dimnames(r2_lm_array) <- list(div_names, radii, NULL)
dimnames(pred_val_lm_array) <- list(div_names, radii, NULL, NULL)
dimnames(coef_array) <- list(div_names, radii, NULL)
coef_df <- melt(coef_array, varnames = c('diversity_type', 'radius', 'iteration'))
r2_lm_df <- melt(r2_lm_array, varnames = c('diversity_type', 'radius', 'iteration'))
pred_val_lm_df <- melt(pred_val_lm_array, varnames = c('diversity_type', 'radius', 'iteration', 'x'))
pred_val_lm_df$x <- newx[pred_val_lm_df$x]

# Predicted values
pred_val_lm_quant <- pred_val_lm_df %>%
  group_by(diversity_type, radius, x) %>%
  summarize(pred_y = quantile(value, probs = 0.5, na.rm = TRUE),
            pred_y_q025 = quantile(value, probs = 0.025, na.rm = TRUE),
            pred_y_q975 = quantile(value, probs = 0.975, na.rm = TRUE),
            pred_y_q25 = quantile(value, probs = 0.25, na.rm = TRUE),
            pred_y_q75 = quantile(value, probs = 0.75, na.rm = TRUE),
            pred_y_mean = mean(value, na.rm = TRUE)) %>%
  arrange(diversity_type, radius, x)

# R2s
r2_lm_quant <-r2_lm_df %>%
  group_by(diversity_type, radius) %>%
  summarize(r2 = quantile(value, probs = 0.5, na.rm = TRUE),
            r2_q025 = quantile(value, probs = 0.025, na.rm = TRUE),
            r2_q975 = quantile(value, probs = 0.975, na.rm = TRUE),
            r2_q25 = quantile(value, probs = 0.25, na.rm = TRUE),
            r2_q75 = quantile(value, probs = 0.75, na.rm = TRUE),
            r2_mean = mean(value, na.rm = TRUE)) %>%
  arrange(diversity_type, radius)

# Slopes or coefficients
coef_quant <- coef_df %>%
  group_by(diversity_type, radius) %>%
  summarize(coef = quantile(value, probs = 0.5, na.rm = TRUE),
            coef_q025 = quantile(value, probs = 0.025, na.rm = TRUE),
            coef_q975 = quantile(value, probs = 0.975, na.rm = TRUE),
            coef_q25 = quantile(value, probs = 0.25, na.rm = TRUE),
            coef_q75 = quantile(value, probs = 0.75, na.rm = TRUE),
            coef_mean = mean(value, na.rm = TRUE)) %>%
  arrange(diversity_type, radius)
  
  
save(pred_val_lm_quant, r2_lm_quant, coef_quant, file = file.path(fp, 'fiafitplotdat_pnw.R'))

###############################

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
  geom_point(aes(x = elevation_sd, y = alpha_diversity), alpha = 0.05) +
  geom_ribbon(data = pred_val_quant %>% filter(diversity_type == 'alpha_diversity'), 
              aes(x = x, ymin = pred_y_min, ymax = pred_y_max), fill = 'blue', alpha = 0.3) +
  geom_line(data = pred_val_quant %>% filter(diversity_type == 'alpha_diversity'), 
            aes(x = x, y = pred_y), color = 'blue', size = 1.5)

ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free') +
  geom_point(aes(x = elevation_sd, y = beta_diversity), alpha = 0.05) +
  geom_ribbon(data = pred_val_quant %>% filter(diversity_type == 'beta_diversity'), 
              aes(x = x, ymin = pred_y_min, ymax = pred_y_max), fill = 'blue', alpha = 0.3) +
  geom_line(data = pred_val_quant %>% filter(diversity_type == 'beta_diversity'), 
            aes(x = x, y = pred_y), color = 'blue', size = 1.5)

ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free') +
  geom_point(aes(x = elevation_sd, y = gamma_diversity), alpha = 0.05) +
  geom_ribbon(data = pred_val_quant %>% filter(diversity_type == 'gamma_diversity'), 
              aes(x = x, ymin = pred_y_min, ymax = pred_y_max), fill = 'blue', alpha = 0.3) +
  geom_line(data = pred_val_quant %>% filter(diversity_type == 'gamma_diversity'), 
            aes(x = x, y = pred_y), color = 'blue', size = 1.5)

### LMs

fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/fia_exploratory_plots/USA'

# Variation explained
r2plot <- ggplot(r2_lm_df, aes(x = radius, y = value, group = interaction(radius, diversity_type), fill = diversity_type)) +
  geom_boxplot() +
  theme(legend.position = 'bottom')
ggsave(file.path(fpfig, 'fia_lm_rsquared.png'), r2plot, height = 5, width = 6, dpi = 200)

# Predicted values
pred_val_lm_quant <- pred_val_lm_df %>%
  group_by(diversity_type, radius, x) %>%
  summarize(pred_y = quantile(value, probs = 0.5, na.rm = TRUE),
            pred_y_min = quantile(value, probs = 0.025, na.rm = TRUE),
            pred_y_max = quantile(value, probs = 0.975, na.rm = TRUE)) %>%
  arrange(diversity_type, radius, x)

# Get rid of values out of range.
xranges <- biogeo %>% group_by(radius) %>% summarize_at(dplyr::vars(elevation_sd), funs(sdmin = min, sdmax = max), na.rm = TRUE)
pred_val_lm_quant <- pred_val_lm_quant %>%
  ungroup %>%
  left_join(xranges) %>%
  filter(x <= sdmax)

# Calculate mean (pseudo) r2s for plotting
meanr2s <- r2_lm_df %>%
  group_by(diversity_type, radius) %>%
  summarize(r2 = mean(value),
            r2expr = as.character(eval(substitute(expression(R^2 == x), list(x = round(r2, 2))))))

hexfill <- scale_fill_gradient(low = 'gray90', high = 'black')
hextheme <- theme(strip.background = element_blank())
border <- panel_border(colour = 'black')
radlabel <- labeller(radius = function(x) paste(as.integer(x), 'km'))

alphaplot <- ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = alpha_diversity)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'alpha_diversity'), 
              aes(x = x, ymin = pred_y_min, ymax = pred_y_max), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'alpha_diversity'), 
            aes(x = x, y = pred_y), color = 'red', size = 1.5) +
  geom_text(data = subset(meanr2s, diversity_type == 'alpha_diversity'),
            aes(x = -Inf, y = Inf, label = r2expr),
            parse = TRUE, hjust = -0.5, vjust = 2) +
  hexfill + hextheme + border + labs(x = 'Elevation standard deviation', y = 'Alpha-diversity')

betaplot <- ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = beta_diversity)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'beta_diversity'), 
              aes(x = x, ymin = pred_y_min, ymax = pred_y_max), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'beta_diversity'), 
            aes(x = x, y = pred_y), color = 'red', size = 1.5) +
  geom_text(data = subset(meanr2s, diversity_type == 'beta_diversity'),
            aes(x = Inf, y = -Inf, label = r2expr),
            parse = TRUE, hjust = 1.1, vjust = -2) +
  hexfill + hextheme + border + labs(x = 'Elevation standard deviation', y = 'Beta-diversity')

gammaplot <- ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = gamma_diversity)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'gamma_diversity'), 
              aes(x = x, ymin = pred_y_min, ymax = pred_y_max), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'gamma_diversity'), 
            aes(x = x, y = pred_y), color = 'red', size = 1.5) +
  geom_text(data = subset(meanr2s, diversity_type == 'gamma_diversity'),
            aes(x = -Inf, y = Inf, label = r2expr),
            parse = TRUE, hjust = -0.5, vjust = 2) +
  hexfill + hextheme + border + labs(x = 'Elevation standard deviation', y = 'Gamma-diversity')

ggsave(file.path(fpfig, 'fia_alpha_regressions.png'), alphaplot, height = 4, width = 12, dpi = 300)
ggsave(file.path(fpfig, 'fia_beta_regressions.png'), betaplot, height = 4, width = 12, dpi = 300)
ggsave(file.path(fpfig, 'fia_gamma_regressions.png'), gammaplot, height = 4, width = 12, dpi = 300)
