# Figures for FIA beta-diversity at different scales, regressed on elevational diversity

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'

bd <- read.csv(file.path(fp, 'fia_betatd.csv'), stringsAsFactors = FALSE)
ed <- read.csv(file.path(fp, 'fia_elev_stats_noalaska.csv'), stringsAsFactors = FALSE)

library(dplyr)

bd <- bd %>%
  mutate(radius = radius/1000) %>% # put radius in km
  left_join(ed)

# Scatterplots

library(cowplot)

fia_beta_plot <- ggplot(bd %>% subset(radius %in% c(5,10,20,50)), aes(x = sd_elev, y = beta_pairwise_abundance)) +
  geom_point(alpha = 0.33) +
  stat_smooth() +
  facet_wrap(~ radius, labeller = labeller(radius = function(x) paste(x, 'km')), scales = 'free_x') +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Pairwise abundance-weighted beta-diversity', expand = c(0,0)) +
  labs(x = 'Standard deviation of elevation')

fpfig <- 'C:/Users/Q/Google Drive/NASABiodiversityWG/Figures/betadiv/'
ggsave(file.path(fpfig, 'fiabetadiversity.png'), fia_beta_plot, height = 9, width = 9, dpi = 400)

# Alpha diversity (taxonomic)

load(file.path(fp, 'fia_diversitymetrics.RData'))

# Calculate plot diversity within radii

library(sp)

fiapnw <- read.csv(file.path(fp, 'finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = F)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

plot_diversity <- left_join(plot_diversity, fiacoords)
radii <- c(5, 10, 20, 50)


neighbordiv <- function(x) {
  neighbordists <- spDistsN1(pts = cbind(plot_diversity$lon, plot_diversity$lat), pt = c(x$lon, x$lat), longlat = TRUE)
  commdat <- list()
  for (i in 1:length(radii)) {
    neighbors_incircle <- plot_diversity[neighbordists <= radii[i], ]
    commdat[[i]] <- with(neighbors_incircle, c(radius = radii[i], richness = median(richness,na.rm=T), shannon_basalarea = median(shannon_basalarea, na.rm=T), mpd = median(mpd_z, na.rm=T), mntd = median(mntd_z, na.rm=T)))
  }
  as.data.frame(do.call('rbind', commdat))
}

fia_alpha <- plot_diversity %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>% 
  do(neighbordiv(.)) # Takes 4 minutes on my machine.


# Join with elevational diversity
ad <- left_join(fia_alpha, ed)

fia_alpha_plot <- ggplot(ad, aes(x = sd_elev, y = shannon_basalarea)) +
  geom_point(alpha = 0.33) +
  stat_smooth() +
  facet_wrap(~ radius, labeller = labeller(radius = function(x) paste(x, 'km')), scales = 'free_x') +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Median Shannon alpha-diversity', expand = c(0,0)) +
  labs(x = 'Standard deviation of elevation')

fpfig <- 'C:/Users/Q/Google Drive/NASABiodiversityWG/Figures/betadiv/'
ggsave(file.path(fpfig, 'fiaalphadiversity.png'), fia_alpha_plot, height = 9, width = 9, dpi = 400)
