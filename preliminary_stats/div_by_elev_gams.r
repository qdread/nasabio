# Generalized additive model fit statistics for fia/bbs alpha and beta diversity at different radii.
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
  do(fit = gam(shannon_basalarea ~ sd_elev, data = .))

alpha_fitdat <- alphagams %>% summarize(rsq = summary(fit)$r.sq) %>% cbind(radius = radii)

betagams <- bd %>%
  filter(radius %in% radii) %>%
  group_by(radius) %>%
  do(fit = gam(beta_pairwise_abundance ~ sd_elev, data = .))

beta_fitdat <- betagams %>% summarize(rsq = summary(fit)$r.sq) %>% cbind(radius = radii)

fiafitdat <- rbind(data.frame(diversity='alpha',alpha_fitdat), data.frame(diversity='beta',beta_fitdat))

library(cowplot)
pgamfia <- ggplot(fiafitdat, aes(x=factor(radius), y=rsq)) + 
  geom_point(size = 3) + 
  facet_wrap(~ diversity) +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  labs(x = 'Radius (km)', y = expression(R^2)) +
  ggtitle('FIA: GAM fits by radius')
ggsave(file.path(fpfig, 'gam_fits_by_radius_FIA.png'), pgamfia, height = 4, width = 7, dpi = 400)

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

bbsfitdat <- rbind(data.frame(diversity='alpha',alpha_fitdat), data.frame(diversity='beta',beta_fitdat))

pgambbs <- ggplot(bbsfitdat, aes(x=factor(radius), y=rsq)) + 
  geom_point(size = 3) + 
  facet_wrap(~ diversity) +
  panel_border(colour = 'black') +
  theme(strip.background = element_blank()) +
  labs(x = 'Radius (km)', y = expression(R^2)) +
  ggtitle('BBS: GAM fits by radius')
ggsave(file.path(fpfig, 'gam_fits_by_radius_BBS.png'), pgambbs, height = 4, width = 7, dpi = 400)