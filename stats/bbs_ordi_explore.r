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
# Also use only the geodiversity variables that have 5k resolution for better comparison.
bbsbio <- bbsbio %>%
  select(rteNo, lon, lat, lon_aea, lat_aea, contains('point'), contains('_100'))

bbsgeo <- bbsgeo %>%
  select(rteNo, lat, lon, HUC4, contains('point'), matches('5k.*100'))


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
pred_category[grep('geo|soil', pred_names)] <- 'geology'
# Split up bioclim variables into those about precip and those about temp
# temp 1-11, precip 12-19
for (i in 1:19) {
  x <- ifelse(i <=11, 'temp', 'precip')
  pred_category[grep(paste0('bio',i,'_'), pred_names)] <- x
}

pred_domain <- rep('local', length(pred_names))
pred_domain[grep('mean', pred_names)] <- 'mean'
pred_domain[grep('sd|diversity|richness', pred_names)] <- 'sd'
pred_domain[grepl('tri', pred_names) & !grepl('point', pred_names)] <- 'tri'
pred_domain[grepl('rough', pred_names) & !grepl('point', pred_names)] <- 'roughness'

pred_df <- data.frame(variable = pred_names,
                            domain = pred_domain,
                            category = pred_category)

# Calculate how each axis is weighted on the different categories.
# each column of rotation is the axis loadings
# do one through 10

library(reshape2)
library(ggplot2)

pred_df <- cbind(pred_df, bbsgeopc$rotation[,1:10])

pred_byaxis <- pred_df %>%
  select(-variable) %>%
  melt(id.vars = c('domain','category'), value.name = 'loading') %>%
  group_by(variable, domain, category) %>%
  summarize(loading = sum(abs(loading)))

ggplot(pred_byaxis, aes(x = variable, y = loading)) +
  geom_col(aes(fill = category), position = 'stack') +
  theme_classic()

ggplot(pred_byaxis, aes(x = variable, y = loading)) +
  geom_col(aes(fill = domain), position = 'stack') +
  theme_classic()

# Try to look at the top few variables for each predictor.

pred_melt <- melt(pred_df, id.vars = c('variable','domain','category'), variable.name = 'axis', value.name = 'loading')

get_top <- function(x, ntop = 3) {
  x[order(abs(x$loading), decreasing = TRUE), ][1:ntop, c('variable','domain','category', 'loading')]
}

pred_top <- pred_melt %>%
  group_by(axis) %>%
  do(get_top(., ntop = 5))

# We can get (sort of) axis interpretations from this:
# 1 Heterogeneity of climate in region (temperature and cloud cover)
# 2 A mixed bag climate axis
# 3 Local and mean temperature
# 4 Cloud cover and precip locally
# 5 Local heterogeneity of temperature
# 6 Local heterogeneity of precip
# 7 Human impacts
# 8 Local topography
# 9 Heterogeneity of climate in region
# 10 local heterogeneity of cloud cover

# together these ten axes explain over 70% of the variation in "geodiversity".
# they are really slanted toward climate

# Use geodiversity to predict biodiversity --------------------------------

# Options: use reduced things to predict, or use other techniques like lasso with shrinkage.

# Lasso
library(glmnet)

# Predict alpha (local), beta over 100 km, and gamma over 100 km
bbs3bio <- bbsbio %>%
  select(rteNo, lon, lat, alpha_richness_point, beta_td_sorensen_pa_100, gamma_richness_100) %>%
  rename(alpha = alpha_richness_point, beta = beta_td_sorensen_pa_100, gamma = gamma_richness_100)

# Aspect has a lot of missing values so get rid.

bbsalpha_dat <- left_join(bbs3bio[,c('rteNo','alpha')], bbsgeo[,-(2:4)]) %>% 
  select(-contains('mode'), -contains('aspect')) %>%
  filter(complete.cases(.))

# Select training set
n <- nrow(bbsalpha_dat)
set.seed(2002)
train_idx <- sample(n, n/2, replace = FALSE)

lambdas <- 10 ^ seq(10, -2, length.out = 100)

dat_x <- as.matrix(bbsalpha_dat[train_idx, -(1:2)])
dat_y <- bbsalpha_dat$alpha[train_idx]
test_y <- bbsalpha_dat$alpha[-train_idx]
lasso_mod <- glmnet(x = dat_x, y = dat_y, alpha = 1, lambda = lambdas)

set.seed(2003)
cv_out <- cv.glmnet(x = dat_x, y = dat_y, alpha = 1)

# With full model
set.seed(2004)
cv_out <- cv.glmnet(x = as.matrix(bbsalpha_dat[,-(1:2)]), y = bbsalpha_dat$alpha, alpha = 1)

# Compare to linear model
alpha_glm <- glm(alpha ~ ., data = bbsalpha_dat[,-1], family = 'gaussian')
library(boot) 
cv_alpha_glm <- cv.glm(bbsalpha_dat[,-1], alpha_glm, K = 10) # tenfold cross validation shows that lasso model is better than full.

# Best lambda is chosen to be near zero
lasso_mod <- glmnet(x = as.matrix(bbsalpha_dat[,-(1:2)]), y = bbsalpha_dat$alpha, alpha = 1, lambda = lambdas)
lasso_pred <- predict(lasso_mod, s = cv_out$lambda.min, newx = as.matrix(bbsalpha_dat[, -(1:2)]))
plot(bbsalpha_dat$alpha ~ lasso_pred)
sqrt(mean((lasso_pred - test_y)^2)) #RMSE, it's off by 14.
abline(0,1,col='red')

lasso_coef <- predict(lasso_mod, type = 'coefficients', s = cv_out$lambda.min) # 97 of the 248 coefficients are not zero. If all data are used with 10 fold CV, there are still 183 coefficients preserved.

### 
# Prediction of beta diversity

bbsbeta_dat <- left_join(bbs3bio[,c('rteNo','beta')], bbsgeo[,-(2:4)]) %>% 
  select(-contains('mode'), -contains('aspect')) %>%
  filter(complete.cases(.))

# Select training set
n <- nrow(bbsbeta_dat)
set.seed(5005)
train_idx <- sample(n, ceiling(n/2), replace = FALSE)

lambdas <- 10 ^ seq(10, -2, length.out = 100)

dat_x <- as.matrix(bbsbeta_dat[train_idx, -(1:2)])
dat_y <- qlogis(bbsbeta_dat$beta[train_idx])
test_y <- qlogis(bbsbeta_dat$beta[-train_idx])
lasso_mod <- glmnet(x = dat_x, y = dat_y, alpha = 1, lambda = lambdas)

set.seed(5006)
cv_out <- cv.glmnet(x = dat_x, y = dat_y, alpha = 1)

# Best lambda is chosen to be near zero
lasso_pred <- predict(lasso_mod, s = cv_out$lambda.min, newx = as.matrix(bbsbeta_dat[-train_idx, -(1:2)]))
plot(plogis(test_y) ~ plogis(lasso_pred)) # Back transformation
sqrt(mean((plogis(lasso_pred) - plogis(test_y))^2)) # RMSE
abline(0,1,col='red') 

lasso_coef <- predict(lasso_mod, type = 'coefficients', s = cv_out$lambda.min) # 38 of the 248 coefficients are not zero.


###
# predict gamma
bbsgamma_dat <- left_join(bbs3bio[,c('rteNo','gamma')], bbsgeo[,-(2:4)]) %>% 
  select(-contains('mode'), -contains('aspect')) %>%
  filter(complete.cases(.))

# Select training set
n <- nrow(bbsgamma_dat)
set.seed(707)
train_idx <- sample(n, n/2, replace = FALSE)

lambdas <- 10 ^ seq(10, -2, length.out = 100)

dat_x <- as.matrix(bbsgamma_dat[train_idx, -(1:2)])
dat_y <- bbsgamma_dat$gamma[train_idx]
test_y <- bbsgamma_dat$gamma[-train_idx]
lasso_mod <- glmnet(x = dat_x, y = dat_y, alpha = 1, lambda = lambdas)

set.seed(708)
cv_out <- cv.glmnet(x = dat_x, y = dat_y, alpha = 1)

# Best lambda is a bit bigger.
lasso_pred <- predict(lasso_mod, s = cv_out$lambda.min, newx = as.matrix(bbsalpha_dat[-train_idx, -(1:2)]))
plot(test_y ~ lasso_pred)
sqrt(mean((lasso_pred - test_y)^2)) #RMSE, it's off by 14.
abline(0,1,col='red')

lasso_coef <- predict(lasso_mod, type = 'coefficients', s = cv_out$lambda.min) # 132 of the 248 coefficients are not zero.



####
# Tenfold cross validation.

