# Simulation of radii

# Generate a bunch of random x,y coordinates, each with a value of geo and of bio, with a linear relationship.

b0 <- 10
b1 <- 5

set.seed(555)

# Make an underlying grid of "geodiversity" decreasing as you go away from the origin (representing true elevation SD at that point, no scale dependence)
# Then generate random true geo values for points on the grid chosen from a uniform distribution
# The geodiversity value has the same mean but SD is dependent on the distance from the origin
# Biodiversity is linearly related to that geodiversity

sd_max <- 10

# sd gets smaller as you go away from origin
sd_distance <- function(d, sd0) sd0 - .01 * d

# Create tiled surface of values with that sd
geo_mat <- matrix(NA, nrow=1000, ncol=1000)
sd_mat <- matrix(NA, nrow=1000, ncol=1000)
for (i in 1:1000) {
  for (j in 1:1000) {
    sd_mat[i,j] <- sd_distance(sqrt((i-500)^2 + (j-500)^2), sd_max)
    geo_mat[i,j] <- rnorm(1, mean = 0, sd = sd_mat[i,j])
  }
}

# plot to show
library(raster)
geo_raster <- raster(geo_mat)
plot(geo_raster)
plot(raster(sd_mat))

set.seed(666)

x <- runif(1000, 0, 1000)
y <- runif(1000, 0, 1000)
bio <- numeric(1000)
for (i in 1:1000) {
  bio[i] <- b0 + b1 * sd_mat[ceiling(x[i]), ceiling(y[i])] + rnorm(1000, 0, 10)
}
dat <- data.frame(x, y, bio)

# Calculate the neighborhood averages of bio and geo

radii <- c(5, 10, 20, 50, 100)

bio_avg <- matrix(NA, nrow = nrow(dat), ncol = length(radii))
geo_avg <- matrix(NA, nrow = nrow(dat), ncol = length(radii))
n_points <- matrix(NA, nrow = nrow(dat), ncol = length(radii))

biggrid <- expand.grid(x=1:1000,y=1:1000)

for (i in 1:nrow(dat)) {
  print(i)
  dist_i <- sqrt((dat$x[i] - dat$x)^2 + (dat$y[i] - dat$y)^2)
  grid_dist_i <- sqrt((dat$x[i] - biggrid$x)^2 + (dat$y[i] - biggrid$y)^2)
  for (j in 1:length(radii)) {
    idx_ij <- which(dist_i <= radii[j])
    bio_avg[i, j] <- mean(dat$bio[idx_ij])
    # Calculate geodiversity inside the radius. (sd of all the geo pixel values within the radius)
    grididx_ij <- which(grid_dist_i <= radii[j])
    geo_pts <- numeric(0)
    for (k in 1:length(grididx_ij)) {
      geo_pts[k] <- geo_mat[biggrid$x[grididx_ij[k]], biggrid$y[grididx_ij[k]]]
    }
    geo_avg[i, j] <- sd(geo_pts)
    n_points[i, j] <- length(idx_ij)
  }
}

# Create plots

lms <- list()

for (i in 1:length(radii)) {
  lms[[i]] <- lm(bio_avg[,i] ~ geo_avg[,i])
}

lmsumm <- lapply(lms, summary)
data.frame(radius = radii,
           slope = sapply(lms, function(x) coef(x)[1]),
           r.squared = sapply(lmsumm, function(x) x$r.squared))

par(mfrow = c(1,5))

for (i in 1:length(radii)) {
  plot(x = geo_avg[,i], y = bio_avg[,i], main = paste('radius', radii[i]), xlab = 'geo st dev', ylab = 'bio mean', ylim = c(0, 80), xlim = c(2, 10))
  abline(a = coef(lms[[i]])[1], b = coef(lms[[i]])[2], lwd = 2, col = 'red')
  mtext(paste('r2 = ', round(summary(lms[[i]])$r.sq, 2)), side = 1, line = 4.5)
  
}
