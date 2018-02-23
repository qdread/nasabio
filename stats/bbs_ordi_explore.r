# Exploratory data and visualizations for BBS
# QDR Nasabioxgeo 22 Feb 2018


# Load data ---------------------------------------------------------------

library(dplyr)

# Change file path depending on whether run locally or remotely
#fp <- '/mnt/research/nasabio/data/bbs'
fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'

bbsbio <- read.csv(file.path(fp, 'bbs_allbio_wide.csv'), stringsAsFactors = FALSE)
bbsgeo <- read.csv(file.path(fp, 'bbs_allgeo_wide.csv'), stringsAsFactors = FALSE)

# Code the discrete variables as such.
disc_vars <- grep('mode', names(bbsgeo), value = TRUE)

bbsgeo <- bbsgeo %>%
  mutate_at(disc_vars, as.factor)

# Get rid of all columns except the point values and the 100 km radius values.
bbsbio <- bbsbio %>%
  select(rteNo, lon, lat, lon_aea, lat_aea, contains('point'), contains('_100'))

bbsgeo <- bbsgeo %>%
  select(rteNo, lat, lon, HUC4, contains('point'), contains('_100_'))

# Ordinations on bbs geodiversity -----------------------------------------

# Correlation matrix including the factors.
#library(polycor)
#bbsgeocor <- hetcor(bbsgeo[,-(1:4)])

# For now ignore the factors as there are only 3

bbsgeocor <- bbsgeo %>%
  select(-rteNo, -lat, -lon, -HUC4, -contains('mode')) %>%
  cor(use = 'pairwise.complete.obs')

good_rtes <- bbsgeo$rteNo[complete.cases(bbsgeo)]

# Do principal components 
bbsgeopc <- bbsgeo %>%
  filter(complete.cases(.)) %>%
  select(-rteNo, -lat, -lon, -HUC4, -contains('mode')) %>%
  prcomp(center = TRUE, scale = TRUE)

summ_pc <- summary(bbsgeopc)
summ_pc$importance[, 1:20]

# Look at what the loadings are on the different axes
# Color the axes by categories and then plot them as some kind of heat map.

pred_names <- dimnames(bbsgeopc$rotation)[[1]]
pred_category <- rep('topography', length(pred_names))
pred_category[grep('dhi', pred_names)] <- 'DHI'
pred_category[grep('cloud', pred_names)] <- 'clouds'
pred_category[grep('human|light', pred_names)] <- 'humans'
# Split up bioclim variables into those about precip and those about temp
# temp 1-11, precip 12-19
for (i in 1:19) {
  x <- ifelse(i <=11, 'temp', 'precip')
  pred_category[grep(paste0('bio',i,'_'), pred_names)] <- x
}

pred_domain <- rep('local', length(pred_names))
pred_domain[grep('mean', pred_names)] <- 'mean'
pred_domain[grep('sd', pred_names)] <- 'sd'
pred_domain[grepl('tri', pred_names) & !grepl('point', pred_names)] <- 'tri'
pred_domain[grepl('rough', pred_names) & !grepl('point', pred_names)] <- 'roughness'

pred_df <- data.frame(variable = pred_names,
                            domain = pred_domain,
                            category = pred_category)

# Calculate how each axis is weighted on the different categories.
# each column of rotation is the axis loadings
# do one through 10
pred_df <- cbind(pred_df, bbsgeopc$rotation[,1:10])

pred_byaxis <- pred_df %>%
  select(-variable) %>%
  melt(id.vars = c('domain','category'), value.name = 'loading') %>%
  group_by(variable, domain, category) %>%
  summarize(loading = sum(abs(loading)))
