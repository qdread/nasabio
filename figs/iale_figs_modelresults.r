# Model result figures for IALE talk
# 06 April 2018

# Load coefficient and r squared files

bbs_coefs <- read.csv('C:/Users/Q/google_drive/NASABiodiversityWG/Conferences/IALE2018/bbs_coefs_iale.csv', stringsAsFactors = FALSE)
fia_coefs <- read.csv('C:/Users/Q/google_drive/NASABiodiversityWG/Conferences/IALE2018/fia_coefs_iale.csv', stringsAsFactors = FALSE)
all_coefs <- cbind(taxon = rep(c('bird','tree'), c(nrow(bbs_coefs), nrow(fia_coefs))),
                   rbind(bbs_coefs, fia_coefs)) %>%
  mutate(predictor = factor(predictor, levels = c('sd', 'roughness', 'tri')))


library(ggplot2)

# Slope of precipitation
th <- theme_bw() + theme(strip.background = element_rect(fill = NA),
                         legend.position = 'bottom',
                         panel.grid.minor.x = element_blank())
zeroline <- geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 1)
taxon_colors <- c('goldenrod', 'forestgreen')
metric_colors <- c('#1b9e77','#d95f02','#7570b3')
cs1 <- scale_color_manual(values = taxon_colors)
cs2 <- scale_color_manual(values = metric_colors, name = 'metric', labels = c('standard deviation','roughness','TRI'))
scale_x <- scale_x_continuous(breaks = c(5, 20, 50, 100, 200), labels = c(5, 20, 50, 100, 200))

# Coefficient plots colored by taxon

ggplot(all_coefs, aes(x = radius, y = slope_precip, color = taxon, group = taxon)) +
  zeroline +
  geom_line() + geom_point() +
  facet_grid(response ~ predictor, scales = 'free_y') +
  labs(x = 'radius (km)', y = 'standardized slope (precipitation)') +
  th + cs1 + scale_x

ggplot(all_coefs, aes(x = radius, y = slope_elev, color = taxon, group = taxon)) +
  zeroline +
  geom_line() + geom_point() +
  facet_grid(response ~ predictor, scales = 'free_y') +
  labs(x = 'radius (km)', y = 'standardized slope (elevation)') +
  th + cs1 + scale_x

ggplot(all_coefs, aes(x = radius, y = slope_footprint, color = taxon, group = taxon)) +
  zeroline +
  geom_line() + geom_point() +
  facet_grid(response ~ predictor, scales = 'free_y') +
  labs(x = 'radius (km)', y = 'standardized slope (human footprint)') +
  th + cs1 + scale_x

# Coefficient plots colored by metric

pcoef_precip <- ggplot(all_coefs, aes(x = radius, y = slope_precip, color = predictor, group = predictor)) +
  zeroline +
  geom_line() + geom_point() +
  facet_grid(response ~ taxon, scales = 'free_y') +
  labs(x = 'radius (km)', y = 'standardized slope (precipitation)') +
  th + cs2 + scale_x

pcoef_elev <- ggplot(all_coefs, aes(x = radius, y = slope_elev, color = predictor, group = predictor)) +
  zeroline +
  geom_errorbar(aes(ymin = slope_elev - 1.96*slope_stderr, ymax = slope_elev + 1.96*slope_stderr), color = 'gray50', width = 1) +
  geom_line() + geom_point() +
  facet_grid(response ~ taxon, scales = 'free_y') +
  labs(x = 'radius (km)', y = 'standardized slope (elevation)') +
  th + cs2 + scale_x

pcoef_hf <- ggplot(all_coefs, aes(x = radius, y = slope_footprint, color = predictor, group = predictor)) +
  zeroline +
  geom_line() + geom_point() +
  facet_grid(response ~ taxon, scales = 'free_y') +
  labs(x = 'radius (km)', y = 'standardized slope (human footprint)') +
  th + cs2 + scale_x

# R-squared plots
ggplot(all_coefs, aes(x = radius, y = r2, color = taxon, group = taxon)) +
  geom_line() + geom_point() +
  facet_grid(response ~ predictor) +
  th + cs1 + scale_x

pr2 <- ggplot(all_coefs, aes(x = radius, y = r2, color = predictor, group = predictor)) +
  geom_line() + geom_point() +
  facet_grid(response ~ taxon) +
  labs(x = 'radius (km)', y = expression(R^2)) +
  th + cs2 + scale_x

ggsave(file.path(fpfig, 'coefplot_precip.png'), pcoef_precip, height = 6.5, width = 5, dpi = 400)
ggsave(file.path(fpfig, 'coefplot_elev.png'), pcoef_elev, height = 6.5, width = 5, dpi = 400)
ggsave(file.path(fpfig, 'coefplot_footprint.png'), pcoef_hf, height = 6.5, width = 5, dpi = 400)
ggsave(file.path(fpfig, 'rsquaredplot.png'), pr2, height = 6.5, width = 5, dpi = 400)


ggsave(file.path(fpfig, 'justelev_coefplot.png'), pcoef_elev, height = 6.5*.8, width = 5*.8, dpi = 400)
ggsave(file.path(fpfig, 'justelev_rsquaredplot.png'), pr2, height = 6.5*.8, width = 5*.8, dpi = 400)
