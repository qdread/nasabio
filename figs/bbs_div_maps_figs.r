# Maps and figures for BBS plots.
# Using most recent biodiversity and geodiversity data
# QDR 14 Sep 2017: NASABIOXGEO project
# Modified 15 Sep 2017: added DHI and night light index.
# Modified 21 Sep 2017: added 75km and 150km radius to the maps.
# Modified 7 Nov 2017: replace with the diversity pooled across years.
# Modified 5 Dec 2017: remake biogeo regressions and maps with corrected data; replace multisite with median pairwise beta-diversity

# Load data ---------------------------------------------------------------

fpdata <- 'C:/Users/Q/Dropbox/projects/nasabiodiv' # Can be changed to /mnt/research/nasabio/data/bbs to get data from HPCC
fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_diversity_maps/' # Directory to save images

bd <- read.csv(file.path(fpdata, 'bbs_betatdpdfd_1year.csv'), stringsAsFactors = FALSE) # older beta-div, with td pd and fd
bdpart <- read.csv(file.path(fpdata, 'bbs_betapart_1year.csv'), stringsAsFactors = FALSE) # newer beta-div, with partitioning but no fd
ad <- read.csv(file.path(fpdata, 'bbs_alpha_1year.csv'), stringsAsFactors = FALSE) # alpha (averaged by radius)
gd <- read.csv(file.path(fpdata, 'bbs_gamma_1year.csv'), stringsAsFactors = FALSE) # gamma
ed <- read.csv(file.path(fpdata, 'bbs_geodiversity_stats.csv'), stringsAsFactors = FALSE) # all geodiversity

# Maps --------------------------------------------------------------------

library(cowplot)
library(dplyr)
source('figs/bbs_map_drawing_fns.r')

radii <- c(50000, 75000, 100000, 150000, 200000, 300000)

library(reshape2)

beta_td <- bdpart %>%
  filter(diversity == 'taxonomic', radius %in% radii, !is.na(beta)) %>%
  select(-diversity) %>%
  dcast(rteNo + lon + lat + radius ~ partition) %>%
  mutate(prop_nested = nestedness/total,
         prop_replace = replacement/total)

img_height <- 10
  
draw_bbs_map(dat = beta_td, zvar = 'total',  
                  colscale = scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.75)(9), breaks = c(0,.5,1), limits=c(0,1)),
                  maptitle = 'BBS beta-diversity, taxonomic, incidence-based, total',
                  fp = fpfig, fname = 'beta_div_tax_total.png', img_h = img_height)

draw_bbs_map(dat = beta_td, zvar = 'prop_replace', 
                  colscale = scale_colour_gradientn(name = 'Proportion due to\nreplacement', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9), breaks = c(0,.5,1), limits=c(0,1)),
                  maptitle = 'BBS beta-diversity, taxonomic, incidence-based, proportion due to species replacement',
                  fp = fpfig, fname = 'beta_div_tax_replacement.png', img_h = img_height)

draw_bbs_map(dat = beta_td, zvar = 'prop_nested',
                  colscale = scale_colour_gradientn(name = 'Proportion due to\nnestedness', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9), breaks = c(0,.5,1), limits=c(0,1)),
                  maptitle = 'BBS beta-diversity, taxonomic, incidence-based, proportion due to species nestedness',
                  fp = fpfig, fname = 'beta_div_tax_nestedness.png', img_h = img_height)

beta_pd <- bdpart %>%
  filter(diversity == 'phylogenetic', radius %in% radii, !is.na(beta)) %>%
  select(-diversity) %>%
  dcast(rteNo + lon + lat + radius ~ partition) %>%
  mutate(prop_nested = nestedness/total,
         prop_replace = replacement/total)

draw_bbs_map(dat = beta_pd, zvar = 'total', 
                  colscale = scale_colour_gradientn(name = 'Phylogenetic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.75)(9), breaks = c(0,.5,1), limits=c(0,1)),
                  maptitle = 'BBS beta-diversity, phylogenetic, incidence-based, total',
                  fp = fpfig, fname = 'beta_div_phy_total.png', img_h = img_height)

draw_bbs_map(dat = beta_pd, zvar = 'prop_replace',  
                  colscale = scale_colour_gradientn(name = 'Proportion due to\nreplacement', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9), breaks = c(0,.5,1), limits=c(0,1)),
                  maptitle = 'BBS beta-diversity, phylogenetic, incidence-based, proportion due to species replacement',
                  fp = fpfig, fname = 'beta_div_phy_replacement.png', img_h = img_height)

draw_bbs_map(dat = beta_pd, zvar = 'prop_nested',  
                  colscale = scale_colour_gradientn(name = 'Proportion due to\nnestedness', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9), breaks = c(0,.5,1), limits=c(0,1)),
                  maptitle = 'BBS beta-diversity, phylogenetic, incidence-based, proportion due to species nestedness',
                  fp = fpfig, fname = 'beta_div_phy_nestedness.png', img_h = img_height)


### "Old" version of beta-diversity

bd %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'beta_pd_pairwise_pa_z',  
               colscale = scale_colour_gradientn(name = 'Phylogenetic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS beta-diversity, phylogenetic, incidence-based (z-score)',
               fp = fpfig, fname = 'beta_div_phy_bettermethod.png', img_h = img_height)

bd %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'beta_td_pairwise_pa',  
               colscale = scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS beta-diversity, taxonomic, incidence-based',
               fp = fpfig, fname = 'beta_div_tax_bettermethod.png', img_h = img_height)

bd %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, beta_fd_pairwise_pa_z > -10) %>%
  draw_bbs_map(zvar = 'beta_fd_pairwise_pa_z',  
               colscale = scale_colour_gradientn(name = 'Functional\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS beta-diversity, functional, incidence-based (z-score; a few outliers excluded)',
               fp = fpfig, fname = 'beta_div_fun.png', img_h = img_height)

### Alpha-diversity maps

# By radius
ad %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
draw_bbs_map(zvar = 'richness',  
             colscale = scale_colour_gradientn(name = 'Taxonomic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
             maptitle = 'BBS alpha-diversity, taxonomic, incidence-based (richness)',
             fp = fpfig, fname = 'alpha_div_tax.png', img_h = img_height)

ad %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'MPD_pa_z',  
               colscale = scale_colour_gradientn(name = 'Phylogenetic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS alpha-diversity, phylogenetic, incidence-based (z-score)',
               fp = fpfig, fname = 'alpha_div_phy.png', img_h = img_height)

ad %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'MPDfunc_pa_z',  
               colscale = scale_colour_gradientn(name = 'Functional\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS alpha-diversity, functional, incidence-based (z-score)',
               fp = fpfig, fname = 'alpha_div_fun.png', img_h = img_height)

### Gamma-diversity maps

# By radius
gd %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'richness',  
               colscale = scale_colour_gradientn(name = 'Taxonomic\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS gamma-diversity, taxonomic, incidence-based (richness)',
               fp = fpfig, fname = 'gamma_div_tax.png', img_h = img_height)

gd %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'MPD_pa_z',  
               colscale = scale_colour_gradientn(name = 'Phylogenetic\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS gamma-diversity, phylogenetic, incidence-based (z-score)',
               fp = fpfig, fname = 'gamma_div_phy.png', img_h = img_height)

gd %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'MPDfunc_pa_z',  
               colscale = scale_colour_gradientn(name = 'Functional\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS gamma-diversity, functional, incidence-based (z-score)',
               fp = fpfig, fname = 'gamma_div_fun.png', img_h = img_height)


### Geodiversity maps

# Standard deviations of the different variables.
# Or richness if that's needed.

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'elevation') %>%
  draw_bbs_map(zvar = 'sd',  
               colscale = scale_colour_gradientn(name = 'Elevation\nstd. dev.', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS geodiversity: elevation',
               fp = fpfig, fname = 'geodiv_elev.png', img_h = img_height)

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'human_footprint') %>%
  draw_bbs_map(zvar = 'sd',  
               colscale = scale_colour_gradientn(name = 'Human footprint\nstd. dev.', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS geodiversity: human footprint index',
               fp = fpfig, fname = 'geodiv_footprint.png', img_h = img_height)

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'nightlight') %>%
  draw_bbs_map(zvar = 'sd',  
               colscale = scale_colour_gradientn(name = 'Night light\nstd. dev.', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS geodiversity: nighttime light intensity',
               fp = fpfig, fname = 'geodiv_nightlight.png', img_h = img_height)

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'dhi_gpp') %>%
  draw_bbs_map(zvar = 'sd',  
               colscale = scale_colour_gradientn(name = 'DHI GPP\nstd. dev.', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS geodiversity: dynamic habitat index GPP',
               fp = fpfig, fname = 'geodiv_dhi_gpp.png', img_h = img_height)

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'dhi_lai8') %>%
  draw_bbs_map(zvar = 'sd',  
               colscale = scale_colour_gradientn(name = 'DHI LAI\nstd. dev.', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS geodiversity: dynamic habitat index LAI',
               fp = fpfig, fname = 'geodiv_dhi_lai.png', img_h = img_height)

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'dhi_ndvi') %>%
  draw_bbs_map(zvar = 'sd',  
               colscale = scale_colour_gradientn(name = 'DHI NDVI\nstd. dev.', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS geodiversity: dynamic habitat index NDVI',
               fp = fpfig, fname = 'geodiv_dhi_ndvi.png', img_h = img_height)

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'geological_age') %>%
  draw_bbs_map(zvar = 'diversity_geodiv',  
               colscale = scale_colour_gradientn(name = 'Geological age\ndiversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.25)(9)),
               maptitle = 'BBS geodiversity: geological age',
               fp = fpfig, fname = 'geodiv_geologicalage.png', img_h = img_height)

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'soil_type') %>%
  draw_bbs_map(zvar = 'diversity_geodiv',  
               colscale = scale_colour_gradientn(name = 'Soil type\ndiversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.5)(9)),
               maptitle = 'BBS geodiversity: soil type',
               fp = fpfig, fname = 'geodiv_soiltype.png', img_h = img_height)

# Bivariate plots ---------------------------------------------------------

radlabel <- labeller(radius = function(x) paste(as.integer(x), 'km'))
radii <- c(50, 75, 100, 150, 200, 300)

# Taxonomic beta-diversity by elevational diversity

td_total_ed_plot <- ed %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  right_join(bd) %>%
  ggplot(aes(x = sd, y = beta_td_pairwise_pa)) +
  geom_point(alpha = 0.05, size = 0.25) +
  stat_smooth(color = 'red', se = FALSE, method = 'auto') +
  facet_grid(. ~ radius, labeller = radlabel) +
  scale_x_continuous(limits=c(0,1250), breaks=c(0,500,1000)) +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Taxonomic beta-diversity', expand = c(0,0), limits = c(0,1)) +
  labs(x = 'Elevation variability')

# td_nested_ed_plot <- ed %>%
#   mutate(radius = radius*1000) %>%
#   filter(variable == 'elevation', radius %in% radii) %>%
#   right_join(beta_td) %>%
#   ggplot(aes(x = sd, y = prop_nested)) +
#   geom_point(alpha = 0.05, size = 0.25) +
#   stat_smooth(color = 'red', se = FALSE, method = 'auto') +
#   facet_grid(. ~ radius, labeller = radlabel) +
#   scale_x_continuous(limits=c(0,1250), breaks=c(0,500,1000)) +
#   theme(strip.background = element_blank()) +
#   panel_border(colour='black') +
#   scale_y_continuous(name = 'Nestedness proportion TD', expand = c(0,0), limits = c(0,1)) +
#   labs(x = 'Elevation variability')

# Phylogenetic beta diversity by elevation diversity

pd_total_ed_plot <- ed %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  right_join(bd) %>%
  ggplot(aes(x = sd, y = beta_pd_pairwise_pa_z)) +
  geom_point(alpha = 0.05, size = 0.25) +
  stat_smooth(color = 'red', se = FALSE, method = 'auto') +
  facet_grid(. ~ radius, labeller = radlabel) +
  scale_x_continuous(limits=c(0,1250), breaks=c(0,500,1000)) +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Phylogenetic beta-diversity (z-score)', expand = c(0,0)) +
  labs(x = 'Elevation variability')

# pd_nested_ed_plot <- ed %>%
#   mutate(radius = radius*1000) %>%
#   filter(variable == 'elevation', radius %in% radii) %>%
#   right_join(beta_pd) %>%
#   ggplot(aes(x = sd, y = prop_nested)) +
#   geom_point(alpha = 0.05, size = 0.25) +
#   stat_smooth(color = 'red', se = FALSE, method = 'auto') +
#   facet_grid(. ~ radius, labeller = radlabel) +
#   scale_x_continuous(limits=c(0,1250), breaks=c(0,500,1000)) +
#   theme(strip.background = element_blank()) +
#   panel_border(colour='black') +
#   scale_y_continuous(name = 'Nestedness proportion PD', expand = c(0,0), limits = c(0,1)) +
#   labs(x = 'Elevation variability')

# Functional beta-diversity by elevation diversity

fd_total_ed_plot <- ed %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  right_join(bd) %>%
  filter(beta_fd_pairwise_pa_z > -10) %>%
  ggplot(aes(x = sd, y = beta_fd_pairwise_pa_z)) +
  geom_point(alpha = 0.05, size = 0.25) +
  stat_smooth(color = 'red', se = FALSE, method = 'auto') +
  facet_grid(. ~ radius, labeller = radlabel) +
  scale_x_continuous(limits=c(0,1250), breaks=c(0,500,1000)) +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Functional beta-diversity (z-score)', expand = c(0,0)) +
  labs(x = 'Elevation variability')

ggsave(file.path(fpfig, 'plot_beta_tax_total_by_elev_sd_bettermethod.png'), td_total_ed_plot, height = 4, width = 12, dpi = 400)
#ggsave(file.path(fpfig, 'plot_beta_tax_nestedness_by_elev_sd.png'), td_nested_ed_plot, height = 4, width = 12, dpi = 400)
ggsave(file.path(fpfig, 'plot_beta_phy_total_by_elev_sd_bettermethod.png'), pd_total_ed_plot, height = 4, width = 12, dpi = 400)
#ggsave(file.path(fpfig, 'plot_beta_phy_nestedness_by_elev_sd.png'), pd_nested_ed_plot, height = 4, width = 12, dpi = 400)
ggsave(file.path(fpfig, 'plot_beta_fun_total_by_elev_sd_bettermethod.png'), fd_total_ed_plot, height = 4, width = 12, dpi = 400)

# Create loop of predictor by response variables and make bivariate plots of all of them.

bbs_xy_plot <- function(xdat, xvar, ydat, yvar, xname, yname, radii, summary_stat = 'sd') {
  xcolumn <- ifelse(xvar %in% c('geological_age', 'soil_type'), 'diversity_geodiv', summary_stat)
  xdat %>%
    filter(variable == xvar) %>%
    right_join(ydat) %>%
    filter(radius %in% radii) %>%
    ggplot(aes_string(x = xcolumn, y = yvar)) +
    geom_hex() +
    scale_fill_gradient(low = 'gray90', high = 'black') +
    stat_smooth(color = 'red', se = FALSE, method = 'auto') +
    facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(x, 'km'))) +
    scale_x_continuous(name = xname) +
    scale_y_continuous(name = yname, expand = c(0,0)) +
    theme(strip.background = element_blank()) +
    panel_border(colour='black')
}

# Loop through all the combinations and make a pdf. This can be used for some crude exploratory data analysis.
# Just use the 1k bioclim.
beta_td <- mutate(beta_td, radius = radius/1000)
beta_pd <- mutate(beta_pd, radius = radius/1000)

# Get rid of some outliers in beta functional diversity
bd$beta_fd_pairwise_pa_z[bd$beta_fd_pairwise_pa_z < -10] <- NA

xvars <- unique(ed$variable)
xvars <- xvars[!grepl('5k', xvars)]
xvars <- xvars[!grepl('bio2|bio3|bio7|bio8|bio9|bio10|bio11|bio16|bio17|bio18|bio19', xvars)]
xvars <- xvars[!grepl('cloud5|cloud6|cloud7|cloud8', xvars)]

ydats <- c('ad', 'ad', 'ad', 'bd', 'bd', 'bd', 'gd', 'gd', 'gd')
yvars <- c('richness', 'MPD_pa_z', 'MPDfunc_pa_z', 'beta_td_pairwise_pa', 'beta_pd_pairwise_pa_z', 'beta_fd_pairwise_pa_z', 'richness', 'MPD_pa_z', 'MPDfunc_pa_z')

xnames <- c('Elevation stdev', 'Slope stdev', 'TPI stdev', 'sin(aspect) stdev', 'cos(aspect) stdev', 'Mean annual temp stdev', 'Temperature seasonality stdev', 'Max temp warmest month stdev', 'Min temp coldest month stdev', 'Mean annual precip stdev', 'Predip wettest month stdev', 'Precip driest month stdev', 'Precip seasonality stdev', 'Mean annual cloudiness stdev', 'Max monthly cloudiness stdev' , 'Min monthly cloudiness stdev', 'Cloudiness seasonality stdev', 'Human footprint index stdev', 'Night light index stdev', 'DHI fPAR stdev', 'DHI GPP stdev', 'DHI LAI stdev', 'DHI NDVI stdev', 'geological age diversity', 'soil type diversity')
ynames <- c('Taxonomic alpha', 'Phylogenetic alpha', 'Functional alpha', 'Taxonomic beta', 'Phylogenetic beta', 'Functional beta', 'Taxonomic gamma', 'Phylogenetic gamma' , 'Functional gamma')

bbs_xy_list <- list()

for (i in 1:length(xvars)) {
  for (j in 1:length(yvars)) {
    bbs_xy_list[[length(bbs_xy_list) + 1]] <- bbs_xy_plot(xdat = ed,
                                                          xvar = xvars[i],
                                                          ydat = get(ydats[j]),
                                                          yvar = yvars[j],
                                                          xname = xnames[i],
                                                          yname = ynames[j],
                                                          radii = c(50, 75, 100, 150, 200, 300))
  }
}

# List of ggplots to multipage pdf
pdf(file.path(fpfig, 'bbs_all_bioxgeo.pdf'), height = 4, width = 12, onefile = TRUE)
  for (i in 1:length(bbs_xy_list)) {
    print(bbs_xy_list[[i]])
  }
dev.off()

### Means of variables (to contrast with standard deviation)

xvars <- xvars[1:23]
xnames <- c('Elevation mean', 'Slope mean', 'TPI mean', 'sin(aspect) mean', 'cos(aspect) mean', 'Mean annual temp mean', 'Temperature seasonality mean', 'Max temp warmest month mean', 'Min temp coldest month mean', 'Mean annual precip mean', 'Predip wettest month mean', 'Precip driest month mean', 'Precip seasonality mean', 'Mean annual cloudiness mean', 'Max monthly cloudiness mean' , 'Min monthly cloudiness mean', 'Cloudiness seasonality mean', 'Human footprint index mean', 'Night light index mean', 'DHI fPAR mean', 'DHI GPP mean', 'DHI LAI mean', 'DHI NDVI mean')

bbs_xymean_list <- list()

for (i in 1:length(xvars)) {
  for (j in 1:length(yvars)) {
    bbs_xymean_list[[length(bbs_xymean_list) + 1]] <- bbs_xy_plot(xdat = ed,
                                                                  xvar = xvars[i],
                                                                  ydat = get(ydats[j]),
                                                                  yvar = yvars[j],
                                                                  xname = xnames[i],
                                                                  yname = ynames[j],
                                                                  radii = c(50, 75, 100, 150, 200, 300),
                                                                  summary_stat = 'mean')
  }
}

# List of ggplots to multipage pdf
pdf(file.path(fpfig, 'bbs_all_bioxgeoMEAN.pdf'), height = 4, width = 12, onefile = TRUE)
for (i in 1:length(bbs_xymean_list)) {
  print(bbs_xymean_list[[i]])
}
dev.off()
