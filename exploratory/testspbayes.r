# Test spBayes with regions

d1 <- bbsbio %>%
  dplyr::select(rteNo, alpha_richness_100) %>%
  left_join(bbsgeo %>% dplyr::select(rteNo, BCR, elevation_5k_100_sd))

bbscoords <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/bbs_correct_route_centroids.csv')

d1 <- left_join(d1, bbscoords[,c(1,4,5)])

d1 <- d1[complete.cases(d1),]

par_start <- list(sigma.sq = 50, phi = 5, tau.sq = 1)
par_priors <- list(sigma.sq.ig = c(2,2), phi.unif = c(1, 100), tau.sq.ig = c(2,0.1))
par_tuning <- list(phi = 0.1, sigma.sq = 0.1, tau.sq = 0.1)

splmtest <- spLM(alpha_richness_100 ~ elevation_5k_100_sd, data = d1, coords = cbind(d1$lon, d1$lat), cov.model = 'gaussian', n.samples = 500, starting = par_start, priors = par_priors, tuning = par_tuning)

splmtest_knots <- spLM(alpha_richness_100 ~ elevation_5k_100_sd, data = d1, coords = cbind(d1$lon, d1$lat), knots = c(5, 5, 1), cov.model = 'gaussian', n.samples = 10000, starting = par_start, priors = par_priors, tuning = par_tuning)

spkpars <- spRecover(splmtest_knots)

spDiag(spkpars)

# Set knots to centroids of the BCR polygons for instance.
library(sp)
library(rgdal)

bcr <- readOGR(dsn = file.path(fp, 'regions'), layer = 'BCR_Terrestrial_master')
bcr <- subset(bcr, COUNTRY %in% 'USA' & !PROVINCE_S %in% 'ALASKA')
library(rgeos)
bnames <- bcr@data$BCRNAME
bcr_combined <- gUnaryUnion(bcr, id = as.character(bnames))

bcr_centroids <- gCentroid(bcr_combined, byid = TRUE)
aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
bcr_cent_aea <- spTransform(bcr_centroids, CRS(aea_crs))

bcr_cents <- as.data.frame(bcr_cent_aea@coords)
bcr_cents$BCR <- dimnames(bcr_cents)[[1]]
names(bcr_cents) <- c('centroid_lon', 'centroid_lat', 'BCR')

splmtest_knots <- spLM(alpha_richness_100 ~ elevation_5k_100_sd, data = d1, coords = cbind(d1$lon, d1$lat), knots = cbind(bcr_cents$centroid_lon, bcr_cents$centroid_lat), cov.model = 'gaussian', n.samples = 5000, starting = par_start, priors = par_priors, tuning = par_tuning)

plot(1:5000, splmtest_knots$p.beta.samples[,1], type = 'l')
plot(1:5000, splmtest_knots$p.beta.samples[,2], type = 'l')

spkpars <- spRecover(splmtest_knots, start = 4001, end = 5000)

spkdiag <- spDiag(spkpars)
