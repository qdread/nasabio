# Map of FIA.


fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'

library(dplyr)
load(file.path(fp, 'fia_diversitymetrics.RData'))

fiapnw <- read.csv(file.path(fp, 'finley_trees_pnw_2015evaluations_feb14_2017.csv'), stringsAsFactors = FALSE)

# Join plot diversity with coordinates
plot_coords <- fiapnw %>% group_by(PLT_CN) %>% summarize(lat = mean(LAT_FUZZSWAP), lon = mean(LON_FUZZSWAP))


plot_diversity <- left_join(plot_diversity, plot_coords)


# Make plots of diversity

library(ggplot2)
latbds <- c(32.5, 50)
lonbds <- c(-125, -114)
aklatbds <- c(52, 61.5)
aklonbds <- c(-153.8, -125)
bothlatbds <- c(32.5, 61.5)
bothlonbds <- c(-153.8, -114)

blackmaptheme <- theme_void() + theme(panel.grid = element_blank(), 
                                      panel.background = element_rect(color = 'black', fill = 'black'), 
                                      plot.background = element_rect(color = 'black', fill = 'black'),  
                                      legend.position = c(0.2, 0.1), 
                                      legend.direction = 'horizontal',
                                      text = element_text(color = 'white'))

purple7colors <- c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')[-(1:2)]
colscale <- scale_colour_gradientn(name = 'Tree species\nrichness', colours = purple7colors)	  
colscalerich <- scale_colour_gradientn(name = 'Tree Shannon\ndiversity', colours = purple7colors)	  


# Continental portion
ggplot(plot_diversity, aes(x = lon, y = lat, color = richness)) +
  borders('state', fill = 'gray70') +
  geom_point() +
  coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = lonbds, ylim = latbds) +
  blackmaptheme + colscale

# Alaska portion
ggplot(plot_diversity, aes(x = lon, y = lat, color = richness)) +
  borders('world', 'canada', fill = 'gray70') +
  borders('world', 'usa', fill = 'gray70') +
  geom_point() +
  coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = aklonbds, ylim = aklatbds) +
  blackmaptheme + colscale

# Both
fiamap1 <- ggplot(plot_diversity, aes(x = lon, y = lat, color = richness)) +
  borders('world', 'canada', fill = 'gray70') +
  borders('world', 'usa', fill = 'gray70') +
  borders('state') +
  geom_point(size = 0.75) +
  coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = bothlonbds, ylim = bothlatbds) +
  blackmaptheme + colscale

# Shannon diversity, basal area
fiamap2 <- ggplot(plot_diversity, aes(x = lon, y = lat, color = shannon_basalarea)) +
  borders('world', 'canada', fill = 'gray70') +
  borders('world', 'usa', fill = 'gray70') +
  borders('state') +
  geom_point(size = 0.75) +
  coord_map(projection = 'albers', lat0=23, lat1=29.5, xlim = bothlonbds, ylim = bothlatbds) +
  blackmaptheme + colscalerich

fp_fig <- 'C:/Users/Q/Google Drive/NasaBiodiversityWG/Figures/descriptivemaps'
ggsave(file.path(fp_fig, 'fia_richness.png'), fiamap1, height = 6, width = 6, dpi = 400)
ggsave(file.path(fp_fig, 'fia_shannon.png'), fiamap2, height = 6, width = 6, dpi = 400)
