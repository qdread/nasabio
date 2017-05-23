# Generalized additive model fit statistics for fia/bbs alpha and beta diversity at different radii.

# Modified 17 May: make axis numbers. Fit a line to it?

# procedure: 
# 1. fit gam to diversity~sd elev at each radius (do this for both alpha and beta)
# 2. get some sort of fit statistic.
# 3. plot it as a function of radius.

# 1. gam fits

library(mgcv)

testgam <- gam(shannon_basalarea ~ sd_elev, data = ad %>% filter(radius == 5))

radii <- c(5,10,20,50,100)

alphagams <- ad %>%
  filter(radius %in% radii) %>%
  group_by(radius) %>%
  do(fit = gam(shannon ~ sd_elev, data = .))

alpha_fitdat <- alphagams %>% summarize(rsq = summary(fit)$r.sq) %>% cbind(radius = radii)

betagams <- bd %>%
  filter(radius %in% radii) %>%
  group_by(radius) %>%
  do(fit = gam(beta_pairwise_abundance ~ sd_elev, data = .))

beta_fitdat <- betagams %>% summarize(rsq = summary(fit)$r.sq) %>% cbind(radius = radii)

gammagams <- gd %>%
  filter(radius %in% radii) %>%
  group_by(radius) %>%
  do(fit = gam(shannon ~ sd_elev, data = .)) 

gamma_fitdat <- gammagams %>% summarize(rsq = summary(fit)$r.sq) %>% cbind(radius = radii)

fiafitdat <- rbind(data.frame(diversity='alpha',alpha_fitdat), 
                   data.frame(diversity='beta',beta_fitdat), 
                   data.frame(diversity='gamma',gamma_fitdat))

# fit a function to each.
lma <- lm(rsq ~ poly(radius,2), data=fiafitdat, subset= diversity=='alpha')
lmb <- lm(rsq ~ poly(radius,2), data=fiafitdat, subset= diversity=='beta')
lmg <- lm(rsq ~ poly(radius,2), data=fiafitdat, subset= diversity=='gamma')

library(cowplot)
pgamfia <- ggplot(fiafitdat, aes(x=radius, y=rsq)) + 
  stat_smooth(method = lm, formula = y ~ x + I(x^2), se=FALSE, size=1, color='red') +
  geom_point(size = 3) + 
  geom_text(data=data.frame(radius=15, rsq=.49, diversity=c('alpha','beta','gamma'), lab=c('R^2 == 0.944','R^2 == .995','R^2 == 0.989')), aes(label=lab), parse=TRUE) +
  facet_wrap(~ diversity) +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  labs(x = 'Radius (km)', y = expression(R^2)) +
  scale_x_continuous(breaks=radii) +
  ggtitle('FIA: GAM fits by radius')
ggsave(file.path(fpfig, 'gam_fits_by_radius_FIA.png'), pgamfia, height = 4, width = 9, dpi = 400)

# do for bbs as well.
radii <- c(50,75,100)

alphagams <- ad %>%
  filter(radius %in% radii) %>%
  group_by(radius) %>%
  do(fit = gam(richness ~ sd_elev, data = .))

alpha_fitdat <- alphagams %>% summarize(rsq = summary(fit)$r.sq) %>% cbind(radius = radii)

betagams <- bd %>%
  filter(radius %in% radii) %>%
  group_by(radius) %>%
  do(fit = gam(beta_td_pairwise_presence ~ sd_elev, data = .))

beta_fitdat <- betagams %>% summarize(rsq = summary(fit)$r.sq) %>% cbind(radius = radii)

gammagams <- gd %>%
  filter(radius %in% radii) %>%
  group_by(radius) %>%
  do(fit = gam(richness ~ sd_elev, data = .))

gamma_fitdat <- gammagams %>% summarize(rsq = summary(fit)$r.sq) %>% cbind(radius = radii)

bbsfitdat <- rbind(data.frame(diversity='alpha',alpha_fitdat), 
                   data.frame(diversity='beta',beta_fitdat), 
                   data.frame(diversity='gamma',gamma_fitdat))

lma <- lm(rsq ~ poly(radius,2), data=bbsfitdat, subset= diversity=='alpha')
lmb <- lm(rsq ~ poly(radius,2), data=bbsfitdat, subset= diversity=='beta')
lmg <- lm(rsq ~ poly(radius,2), data=fiafitdat, subset= diversity=='gamma')


pgambbs <- ggplot(bbsfitdat, aes(x=radius, y=rsq)) + 
  #stat_smooth(method = lm, formula = y ~ x + I(x^2), se=FALSE, size=1, color='red') +
  geom_point(size = 3) + 
  facet_wrap(~ diversity) +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  labs(x = 'Radius (km)', y = expression(R^2)) +
  scale_x_continuous(breaks=radii) +
  ggtitle('BBS: GAM fits by radius')
ggsave(file.path(fpfig, 'gam_fits_by_radius_BBS.png'), pgambbs, height = 4, width = 9, dpi = 400)