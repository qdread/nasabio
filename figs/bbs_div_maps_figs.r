# Maps and figures for BBS plots.
# Using most recent biodiversity and geodiversity data
# QDR 14 Sep 2017: NASABIOXGEO project


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


# Bivariate plots ---------------------------------------------------------


# Maps --------------------------------------------------------------------

library(cowplot)
source('figs/bbs_map_drawing_fns.r')

radii <- c(50000, 100000, 200000, 300000)

library(reshape2)

beta_td <- bdpart_yrmeans %>%
  filter(diversity == 'taxonomic', radius %in% radii, !is.na(beta)) %>%
  select(-diversity) %>%
  dcast(rteNo + lon + lat + radius ~ partition) %>%
  mutate(prop_nested = nestedness/total,
         prop_replace = replacement/total)
  
draw_bbs_map(dat = beta_td, zvar = 'total',  
                  colscale = scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.75)(9), breaks = c(0,.5,1), limits=c(0,1)),
                  maptitle = 'BBS beta-diversity, taxonomic, incidence-based, total',
                  fp = fpfig, fname = 'beta_div_tax_total.png')

draw_bbs_map(dat = beta_td, zvar = 'prop_replace', 
                  colscale = scale_colour_gradientn(name = 'Proportion due to\nreplacement', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9), breaks = c(0,.5,1), limits=c(0,1)),
                  maptitle = 'BBS beta-diversity, taxonomic, incidence-based, proportion due to species replacement',
                  fp = fpfig, fname = 'beta_div_tax_replacement.png')

draw_bbs_map(dat = beta_td, zvar = 'prop_nested',
                  colscale = scale_colour_gradientn(name = 'Proportion due to\nnestedness', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9), breaks = c(0,.5,1), limits=c(0,1)),
                  maptitle = 'BBS beta-diversity, taxonomic, incidence-based, proportion due to species nestedness',
                  fp = fpfig, fname = 'beta_div_tax_nestedness.png')

beta_pd <- bdpart_yrmeans %>%
  filter(diversity == 'phylogenetic', radius %in% radii, !is.na(beta)) %>%
  select(-diversity) %>%
  dcast(rteNo + lon + lat + radius ~ partition) %>%
  mutate(prop_nested = nestedness/total,
         prop_replace = replacement/total)

draw_bbs_map(dat = beta_pd, zvar = 'total', 
                  colscale = scale_colour_gradientn(name = 'Phylogenetic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.75)(9), breaks = c(0,.5,1), limits=c(0,1)),
                  maptitle = 'BBS beta-diversity, phylogenetic, incidence-based, total',
                  fp = fpfig, fname = 'beta_div_phy_total.png')

draw_bbs_map(dat = beta_pd, zvar = 'prop_replace',  
                  colscale = scale_colour_gradientn(name = 'Proportion due to\nreplacement', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9), breaks = c(0,.5,1), limits=c(0,1)),
                  maptitle = 'BBS beta-diversity, phylogenetic, incidence-based, proportion due to species replacement',
                  fp = fpfig, fname = 'beta_div_phy_replacement.png')

draw_bbs_map(dat = beta_pd, zvar = 'prop_nested',  
                  colscale = scale_colour_gradientn(name = 'Proportion due to\nnestedness', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9), breaks = c(0,.5,1), limits=c(0,1)),
                  maptitle = 'BBS beta-diversity, phylogenetic, incidence-based, proportion due to species nestedness',
                  fp = fpfig, fname = 'beta_div_phy_nestedness.png')


### "Old" version of beta-diversity

bdold_yrmeans %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, fd_z > -10) %>%
  draw_bbs_map(zvar = 'fd_z',  
               colscale = scale_colour_gradientn(name = 'Functional\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS alpha-diversity, functional, incidence-based (z-score; a few outliers excluded)',
               fp = fpfig, fname = 'beta_div_fun.png')

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
             fp = fpfig, fname = 'alpha_div_tax.png')

adradius_yrmeans %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'MPD_z',  
               colscale = scale_colour_gradientn(name = 'Phylogenetic\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS alpha-diversity, phylogenetic, incidence-based (z-score)',
               fp = fpfig, fname = 'alpha_div_phy.png')

adradius_yrmeans %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'MPDfunc_z',  
               colscale = scale_colour_gradientn(name = 'Functional\nalpha-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS alpha-diversity, functional, incidence-based (z-score)',
               fp = fpfig, fname = 'alpha_div_fun.png')

### Gamma-diversity maps

# By radius
gd_yrmeans %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'richness',  
               colscale = scale_colour_gradientn(name = 'Taxonomic\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS gamma-diversity, taxonomic, incidence-based (richness)',
               fp = fpfig, fname = 'gamma_div_tax.png')

gd_yrmeans %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'MPD_z',  
               colscale = scale_colour_gradientn(name = 'Phylogenetic\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS gamma-diversity, phylogenetic, incidence-based (z-score)',
               fp = fpfig, fname = 'gamma_div_phy.png')

gd_yrmeans %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii) %>%
  draw_bbs_map(zvar = 'MPDfunc_z',  
               colscale = scale_colour_gradientn(name = 'Functional\ngamma-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS gamma-diversity, functional, incidence-based (z-score)',
               fp = fpfig, fname = 'gamma_div_fun.png')


### Geodiversity maps

# Standard deviations of the different variables.
# Or richness if that's needed.

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'elevation') %>%
  draw_bbs_map(zvar = 'sd',  
               colscale = scale_colour_gradientn(name = 'Elevation\nstd. dev.', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS geodiversity: elevation',
               fp = fpfig, fname = 'geodiv_elev.png')

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'human_footprint') %>%
  draw_bbs_map(zvar = 'sd',  
               colscale = scale_colour_gradientn(name = 'Human footprint\nstd. dev.', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 1)(9)),
               maptitle = 'BBS geodiversity: human footprint index',
               fp = fpfig, fname = 'geodiv_footprint.png')

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'geological_age') %>%
  draw_bbs_map(zvar = 'diversity',  
               colscale = scale_colour_gradientn(name = 'Geological age\ndiversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.25)(9)),
               maptitle = 'BBS geodiversity: geological age',
               fp = fpfig, fname = 'geodiv_geologicalage.png')

ed %>%
  mutate(radius = radius * 1000) %>%
  filter(radius %in% radii, variable == 'soil_type') %>%
  draw_bbs_map(zvar = 'diversity',  
               colscale = scale_colour_gradientn(name = 'Soil type\ndiversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.5)(9)),
               maptitle = 'BBS geodiversity: soil type',
               fp = fpfig, fname = 'geodiv_soiltype.png')
