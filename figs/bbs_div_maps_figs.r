# Maps and figures for BBS plots.
# Using most recent biodiversity and geodiversity data
# QDR 14 Sep 2017: NASABIOXGEO project
# Modified 15 Sep 2017: added DHI and night light index.
# Modified 21 Sep 2017: added 75km and 150km radius to the maps.

# Load data ---------------------------------------------------------------

fpdata <- 'C:/Users/Q/Dropbox/projects/nasabiodiv' # Can be changed to /mnt/research/nasabio/data/bbs to get data from HPCC
fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_diversity_maps/' # Directory to save images

bd <- read.csv(file.path(fpdata, 'bbs_beta_byroute.csv'), stringsAsFactors = FALSE) # older beta-div, with td pd and fd
bdpart <- read.csv(file.path(fpdata, 'bbs_betapart_byroute.csv'), stringsAsFactors = FALSE) # newer beta-div, with partitioning but no fd
ad <- read.csv(file.path(fpdata, 'bbs_alphadiv.csv'), stringsAsFactors = FALSE) # alpha (individual plots)
adradius <- read.csv(file.path(fpdata, 'bbs_alpha.csv'), stringsAsFactors = FALSE) # alpha (averaged by radius)
gd <- read.csv(file.path(fpdata, 'bbs_gammadiv.csv'), stringsAsFactors = FALSE) # gamma
ed <- read.csv(file.path(fpdata, 'bbs_geodiversity_stats.csv'), stringsAsFactors = FALSE) # all geodiversity

# Calculate means of all routes across years 1997-2016

library(dplyr)

bdpart_yrmeans <- bdpart %>%
  filter(year >= 1997) %>%
  group_by(rteNo, lon, lat, radius, diversity, partition) %>%
  summarize(beta = mean(na.omit(beta)))

ad_yrmeans <- ad %>%
  filter(year >= 1997) %>%
  group_by(rteNo, lon, lat) %>%
  summarize(richness = mean(na.omit(richness)),
            MPD_z = mean(na.omit(MPD_pa_z)),
            MNTD_z = mean(na.omit(MNTD_pa_z)),
            MPDfunc_z = mean(na.omit(MPDfunc_pa_z)),
            MNTDfunc_z = mean(na.omit(MNTDfunc_pa_z)))

adradius_yrmeans <- adradius %>%
  filter(year >= 1997) %>%
  group_by(rteNo, lon, lat, radius) %>%
  summarize(richness = mean(na.omit(richness)),
            MPD_z = mean(na.omit(MPD_pa_z)),
            MNTD_z = mean(na.omit(MNTD_pa_z)),
            MPDfunc_z = mean(na.omit(MPDfunc_pa_z)),
            MNTDfunc_z = mean(na.omit(MNTDfunc_pa_z)))

gd_yrmeans <- gd %>%
  filter(year >= 1997) %>%
  group_by(rteNo, lon, lat, radius) %>%
  summarize(richness = mean(na.omit(richness)),
            MPD_z = mean(na.omit(MPD_pa_z)),
            MNTD_z = mean(na.omit(MNTD_pa_z)),
            MPDfunc_z = mean(na.omit(MPDfunc_pa_z)),
            MNTDfunc_z = mean(na.omit(MNTDfunc_pa_z)))

bdold_yrmeans <- bd %>%
  filter(year >= 1997) %>%
  group_by(rteNo, lon, lat, radius) %>%
  summarize(td = mean(na.omit(beta_td_pairwise_pa)),
            pd = mean(na.omit(beta_pd_pairwise_pa)),
            pd_z = mean(na.omit(beta_pd_pairwise_pa_z)),
            fd = mean(na.omit(beta_fd_pairwise_pa)),
            fd_z = mean(na.omit(beta_fd_pairwise_pa_z)))


# Maps --------------------------------------------------------------------

library(cowplot)
source('figs/bbs_map_drawing_fns.r')

radii <- c(50000, 75000, 100000, 150000, 200000, 300000)

library(reshape2)

beta_td <- bdpart_yrmeans %>%
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

beta_pd <- bdpart_yrmeans %>%
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

bdold_yrmeans %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, fd_z > -10) %>%
  draw_bbs_map(zvar = 'fd_z',  
               colscale = scale_colour_gradientn(name = 'Functional\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS alpha-diversity, functional, incidence-based (z-score; a few outliers excluded)',
               fp = fpfig, fname = 'beta_div_fun.png', img_h = img_height)

### Alpha-diversity maps

# Single maps
draw_bbs_map(dat = ad, zvar = 'richness', by_rad = FALSE,
                         colscale = scale_colour_gradientn(name = 'Taxonomic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
                         maptitle = 'BBS alpha-diversity, taxonomic, incidence-based (richness)',
                         fp = fpfig, fname = 'alpha_div_tax.png', write_to_file = TRUE, img_h = 5, img_w = 8.5)

draw_bbs_map(dat = ad, zvar = 'MPD_z', by_rad = FALSE,
             colscale = scale_colour_gradientn(name = 'Phylogenetic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
             maptitle = 'BBS alpha-diversity, phylogenetic, incidence-based (z-score)',
             fp = fpfig, fname = 'alpha_div_phy.png', write_to_file = TRUE, img_h = 5, img_w = 8.5)

# By radius
adradius_yrmeans %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
draw_bbs_map(zvar = 'richness',  
             colscale = scale_colour_gradientn(name = 'Taxonomic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
             maptitle = 'BBS alpha-diversity, taxonomic, incidence-based (richness)',
             fp = fpfig, fname = 'alpha_div_tax.png', img_h = img_height)

adradius_yrmeans %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'MPD_z',  
               colscale = scale_colour_gradientn(name = 'Phylogenetic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS alpha-diversity, phylogenetic, incidence-based (z-score)',
               fp = fpfig, fname = 'alpha_div_phy.png', img_h = img_height)

adradius_yrmeans %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'MPDfunc_z',  
               colscale = scale_colour_gradientn(name = 'Functional\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS alpha-diversity, functional, incidence-based (z-score)',
               fp = fpfig, fname = 'alpha_div_fun.png', img_h = img_height)

### Gamma-diversity maps

# By radius
gd_yrmeans %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'richness',  
               colscale = scale_colour_gradientn(name = 'Taxonomic\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS gamma-diversity, taxonomic, incidence-based (richness)',
               fp = fpfig, fname = 'gamma_div_tax.png', img_h = img_height)

gd_yrmeans %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'MPD_z',  
               colscale = scale_colour_gradientn(name = 'Phylogenetic\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS gamma-diversity, phylogenetic, incidence-based (z-score)',
               fp = fpfig, fname = 'gamma_div_phy.png', img_h = img_height)

gd_yrmeans %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'MPDfunc_z',  
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
  draw_bbs_map(zvar = 'diversity',  
               colscale = scale_colour_gradientn(name = 'Geological age\ndiversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.25)(9)),
               maptitle = 'BBS geodiversity: geological age',
               fp = fpfig, fname = 'geodiv_geologicalage.png', img_h = img_height)

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'soil_type') %>%
  draw_bbs_map(zvar = 'diversity',  
               colscale = scale_colour_gradientn(name = 'Soil type\ndiversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.5)(9)),
               maptitle = 'BBS geodiversity: soil type',
               fp = fpfig, fname = 'geodiv_soiltype.png', img_h = img_height)

# Bivariate plots ---------------------------------------------------------

radlabel <- labeller(radius = function(x) paste(as.integer(x)/1000, 'km'))

# Taxonomic beta-diversity by elevational diversity

td_total_ed_plot <- ed %>%
  mutate(radius = radius*1000) %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  right_join(beta_td) %>%
  ggplot(aes(x = sd, y = total)) +
  geom_point(alpha = 0.05, size = 0.25) +
  stat_smooth(color = 'red', se = FALSE, method = 'auto') +
  facet_grid(. ~ radius, labeller = radlabel) +
  scale_x_continuous(limits=c(0,1250), breaks=c(0,500,1000)) +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Taxonomic beta-diversity', expand = c(0,0), limits = c(0,1)) +
  labs(x = 'Elevation variability')

td_nested_ed_plot <- ed %>%
  mutate(radius = radius*1000) %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  right_join(beta_td) %>%
  ggplot(aes(x = sd, y = prop_nested)) +
  geom_point(alpha = 0.05, size = 0.25) +
  stat_smooth(color = 'red', se = FALSE, method = 'auto') +
  facet_grid(. ~ radius, labeller = radlabel) +
  scale_x_continuous(limits=c(0,1250), breaks=c(0,500,1000)) +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Nestedness proportion TD', expand = c(0,0), limits = c(0,1)) +
  labs(x = 'Elevation variability')

# Phylogenetic beta diversity by elevation diversity

pd_total_ed_plot <- ed %>%
  mutate(radius = radius*1000) %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  right_join(beta_pd) %>%
  ggplot(aes(x = sd, y = total)) +
  geom_point(alpha = 0.05, size = 0.25) +
  stat_smooth(color = 'red', se = FALSE, method = 'auto') +
  facet_grid(. ~ radius, labeller = radlabel) +
  scale_x_continuous(limits=c(0,1250), breaks=c(0,500,1000)) +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Phylogenetic beta-diversity', expand = c(0,0), limits = c(0,1)) +
  labs(x = 'Elevation variability')

pd_nested_ed_plot <- ed %>%
  mutate(radius = radius*1000) %>%
  filter(variable == 'elevation', radius %in% radii) %>%
  right_join(beta_pd) %>%
  ggplot(aes(x = sd, y = prop_nested)) +
  geom_point(alpha = 0.05, size = 0.25) +
  stat_smooth(color = 'red', se = FALSE, method = 'auto') +
  facet_grid(. ~ radius, labeller = radlabel) +
  scale_x_continuous(limits=c(0,1250), breaks=c(0,500,1000)) +
  theme(strip.background = element_blank()) +
  panel_border(colour='black') +
  scale_y_continuous(name = 'Nestedness proportion PD', expand = c(0,0), limits = c(0,1)) +
  labs(x = 'Elevation variability')

ggsave(file.path(fpfig, 'plot_beta_tax_total_by_elev_sd.png'), td_total_ed_plot, height = 4, width = 12, dpi = 400)
ggsave(file.path(fpfig, 'plot_beta_tax_nestedness_by_elev_sd.png'), td_nested_ed_plot, height = 4, width = 12, dpi = 400)
ggsave(file.path(fpfig, 'plot_beta_phy_total_by_elev_sd.png'), pd_total_ed_plot, height = 4, width = 12, dpi = 400)
ggsave(file.path(fpfig, 'plot_beta_phy_nestedness_by_elev_sd.png'), pd_nested_ed_plot, height = 4, width = 12, dpi = 400)

# Create loop of predictor by response variables and make bivariate plots of all of them.

bbs_xy_plot <- function(xdat, xvar, ydat, yvar, xlims, xbreaks, ylims, ybreaks, xname, yname) {
  xdat %>%
    mutate(radius = radius*1000) %>%
    filter(variable == xvar, radius %in% radii) %>%
    right_join(ydat) %>%
    ggplot(aes_string(x = xvar, y = yvar)) +
    geom_point(alpha = 0.05, size = 0.25) +
    stat_smooth(color = 'red', se = FALSE, method = 'auto') +
    facet_grid(. ~ radius, labeller = labeller(radius = function(x) paste(as.integer(x)/1000, 'km'))) +
    scale_x_continuous(name = xname, limits=xlims, breaks=xbreaks) +
    scale_y_continuous(name = yname, limits=ylims, breaks=ybreaks, expand = c(0,0)) +
    theme(strip.background = element_blank()) +
    panel_border(colour='black')
}

