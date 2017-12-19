# Subsampling of BBS points, and model fitting with the subsamples
# Up to 200 km

# NOTE: All needed CSVs are on the hpcc but I have downloaded them locally because it is faster. 
# Change file path to the second file path to get files from hpcc.

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'
#fp <- '/mnt/research/nasabio/data/bbs'

# Load the elevational diversity and abg diversity data
ed <- read.csv(file.path(fp, 'bbs_geodiversity_stats.csv'))
ad <- read.csv(file.path(fp, 'bbs_alpha_1year.csv'))
bd <- read.csv(file.path(fp, 'bbs_betatdpdfd_1year.csv'))
gd <- read.csv(file.path(fp, 'bbs_gamma_1year.csv'))

radii <- c(50, 75, 100, 150, 200)
div_names <- c('alpha_richness','beta_richness','gamma_richness')

library(dplyr)
library(sp)
library(mgcv)

# Function for iterative search.
source('~/GitHub/nasabio/stats/SRS_iterative.r')
# Functions for flagging edge plots
source('~/GitHub/nasabio/stats/spatial_fns.r')

# Combine into a single data frame.
# For now just use elevation SD as the response variable.
biogeo <- ed %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  dplyr::select(rteNo, lon, lat, lon_aea, lat_aea, radius, sd) %>%
  rename(elevation_sd = sd) %>%
  left_join(ad %>% 
              dplyr::select(rteNo, lon, lat, lon_aea, lat_aea, radius, richness) %>% 
              rename(alpha_richness = richness) %>%
              filter(radius %in% radii)) %>%
  left_join(bd %>% 
              dplyr::select(rteNo, lon, lat, lon_aea, lat_aea, radius, beta_td_pairwise_pa) %>% 
              rename(beta_richness = beta_td_pairwise_pa) %>%
              filter(radius %in% radii)) %>%
  left_join(gd %>% 
              dplyr::select(rteNo, lon, lat, lon_aea, lat_aea, radius, richness) %>% 
              rename(gamma_richness = richness) %>%
              filter(radius %in% radii))

# Add latitude and longitudes of route centroids

bbscoords <- read.csv(file.path(fp, 'bbs_correct_route_centroids.csv')) %>% rename(lon_aea=lon.1, lat_aea=lat.1)
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
bbs_aea <- SpatialPoints(coords = bbscoords[,c('lon_aea','lat_aea')], proj4string = CRS(aea_crs))

# Throw out the ones that are within 200 km of Mexico or Canada

bbs_coast200 <- flag_coast_plots(focal_points = bbs_aea, radius = 200e3, border_countries = c('Canada', 'Mexico'))

bbscoords <- cbind(bbscoords, bbs_coast200)

biogeo <- biogeo %>%
  left_join(bbscoords %>% dplyr::select(rteNo, is_edge)) %>%
  filter(!is_edge, complete.cases(.))

length(unique(biogeo$rteNo)) # approx 2700 are not within 200 km of Mexico or Canada.

bbs_aea_noedge <- bbs_aea[!bbscoords$is_edge & bbscoords$rteNo %in% biogeo$rteNo]
bbscoords_noedge <- subset(bbscoords, !is_edge & rteNo %in% biogeo$rteNo)


# Calculate pairwise distances. (only takes ~2 seconds with brute force method)
bbs_dist <- spDists(bbs_aea_noedge, longlat = FALSE)

# Subsampling with linear models ---------------------------------------------------

n_iter <- 999 # Increase later.
subsample_size <- 35
n_max <- 50

library(betareg)

# x-values for predicted y-values
xrange <- range(biogeo$elevation_sd, na.rm=TRUE)
newx <- round(seq(xrange[1], xrange[2], length.out = 50))

r2_lm_array <- array(NA, dim = c(length(div_names), length(radii), n_iter))
pred_val_lm_array <- array(NA, dim = c(length(div_names), length(radii), n_iter, length(newx)))

pb <- txtProgressBar(0, prod(dim(r2_lm_array)), style = 3)
ii <- 0
set.seed(919)

for (k in 1:n_iter) {
  
  # Subsample using the largest radius (200 km)
  sample_idx <- SRS_iterative(focal_points = bbs_aea_noedge, dist_mat = bbs_dist, radius = 200e3, n = n_max)
  
  # If too many are subsampled, take a random subset of the subsample.
  if (length(sample_idx) > subsample_size) {
    plots_in_sample <- bbscoords_noedge$rteNo[sample(sample_idx, subsample_size)]
  } else {
    plots_in_sample <- bbscoords_noedge$rteNo[sample_idx]
  }
  
  for (j in 1:length(radii)) {
    
    dat <- biogeo %>%
      filter(radius == radii[j], rteNo %in% plots_in_sample)
    
    for (i in 1:length(div_names)) {
      ii <- ii + 1
      setTxtProgressBar(pb, ii)

      # Fit lm if alpha or gamma, fit beta-regression (oddly enough) if beta
      if (div_names[i] == 'beta_richness') {
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


# Convert array to data frame
library(reshape2)

dimnames(r2_lm_array) <- list(div_names, radii, NULL)
dimnames(pred_val_lm_array) <- list(div_names, radii, NULL, NULL)
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


# Plot
library(cowplot)

# Plot the r2 values.
r2_lm_quant$radius_plot <- r2_lm_quant$radius + rep(c(-1,0,1), each = 5)

r2_plot <- ggplot(r2_lm_quant, aes(x = radius_plot, group = interaction(radius, diversity_type), color = diversity_type)) +
  geom_segment(aes(xend = radius_plot, y = r2_q25, yend = r2_q75), size = 1) +
  geom_line(aes(y = r2, group = diversity_type), size = 0.5) +
  geom_point(aes(y = r2), size = 2) +
  scale_color_discrete(name = 'Diversity', labels = c('alpha','beta','gamma')) +
  labs(x = 'Radius (km)', y = expression(R^2)) +
  theme(legend.position = c(0.1, 0.9))

# Get rid of values out of range.
xranges <- biogeo %>% group_by(radius) %>% summarize_at(vars(elevation_sd), funs(sdmin = min, sdmax = max), na.rm = TRUE)
pred_val_lm_quant <- pred_val_lm_quant %>%
  ungroup %>%
  left_join(xranges) %>%
  filter(x <= sdmax)

# Median (pseudo) r2s for plotting
r2_lm_quant <- r2_lm_quant %>%
  ungroup %>%
  group_by(diversity_type, radius) %>%
  mutate(r2expr = as.character(eval(substitute(expression(R^2 == x), list(x = round(r2, 2))))))

hexfill <- scale_fill_gradient(low = 'gray90', high = 'black')
hextheme <- theme(strip.background = element_blank())
border <- panel_border(colour = 'black')
radlabel <- labeller(radius = function(x) paste(as.integer(x), 'km'))

alphaplot <- ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = alpha_richness)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'alpha_richness'), 
              aes(x = x, ymin = pred_y_q025, ymax = pred_y_q975), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'alpha_richness'), 
            aes(x = x, y = pred_y), color = 'red', size = 1.5) +
  geom_text(data = subset(r2_lm_quant, diversity_type == 'alpha_richness'),
            aes(x = -Inf, y = Inf, label = r2expr),
            parse = TRUE, hjust = -0.5, vjust = 2) +
  hexfill + hextheme + border + labs(x = 'Elevation standard deviation', y = 'Alpha-diversity')

betaplot <- ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = beta_richness)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'beta_richness'), 
              aes(x = x, ymin = pred_y_q025, ymax = pred_y_q975), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'beta_richness'), 
            aes(x = x, y = pred_y), color = 'red', size = 1.5) +
  geom_text(data = subset(r2_lm_quant, diversity_type == 'beta_richness'),
            aes(x = Inf, y = -Inf, label = r2expr),
            parse = TRUE, hjust = 1.1, vjust = -2) +
  hexfill + hextheme + border + labs(x = 'Elevation standard deviation', y = 'Beta-diversity')

gammaplot <- ggplot(biogeo) + 
  facet_grid(~ radius, scales = 'free', labeller = radlabel) +
  geom_hex(aes(x = elevation_sd, y = gamma_richness)) +
  geom_ribbon(data = pred_val_lm_quant %>% filter(diversity_type == 'gamma_richness'), 
              aes(x = x, ymin = pred_y_q025, ymax = pred_y_q975), fill = 'red', alpha = 0.3) +
  geom_line(data = pred_val_lm_quant %>% filter(diversity_type == 'gamma_richness'), 
            aes(x = x, y = pred_y), color = 'red', size = 1.5) +
  geom_text(data = subset(r2_lm_quant, diversity_type == 'gamma_richness'),
            aes(x = -Inf, y = Inf, label = r2expr),
            parse = TRUE, hjust = -0.5, vjust = 2) +
  hexfill + hextheme + border + labs(x = 'Elevation standard deviation', y = 'Gamma-diversity')

fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_diversity_maps'

ggsave(file.path(fpfig, 'bbs_r2s.png'), r2_plot, height = 5, width = 6, dpi = 300)
ggsave(file.path(fpfig, 'bbs_alpha_regressions.png'), alphaplot, height = 4, width = 12, dpi = 300)
ggsave(file.path(fpfig, 'bbs_beta_regressions.png'), betaplot, height = 4, width = 12, dpi = 300)
ggsave(file.path(fpfig, 'bbs_gamma_regressions.png'), gammaplot, height = 4, width = 12, dpi = 300)