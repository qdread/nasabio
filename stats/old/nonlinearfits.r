dat <- subset(biogeo, radius == 100)[,c('elevation_sd','gamma_diversity')]
gam1 <- gam(I(exp(gamma_diversity)) ~ elevation_sd, data=dat)
lm1 <- lm(I(exp(gamma_diversity)) ~ elevation_sd, data=dat)
gam2 <- gam(I(exp(gamma_diversity)) ~ s(elevation_sd, k=4), data=dat)
gam2_pred <- predict(gam2, data.frame(elevation_sd = seq(0, 1200, 50)))

with(dat, plot(elevation_sd, exp(gamma_diversity)))
lines(x = seq(0,1200,50), y = gam2_pred, type = 'l', lwd = 2, col = 'red')

mod3 <- lm(I(exp(gamma_diversity)) ~ bs(elevation_sd, knots = c(200,700,1200)), data= dat)
mod3_pred <- predict(mod3, data.frame(elevation_sd = seq(0, 1200, 50)))
lines(x = seq(0,1200,50), y = mod3_pred, type = 'l', lwd = 2, col = 'blue')
