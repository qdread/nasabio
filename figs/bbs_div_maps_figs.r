# Maps and figures for BBS plots.
# Using most recent biodiversity and geodiversity data
# QDR 14 Sep 2017: NASABIOXGEO project


# Load data ---------------------------------------------------------------

fpdata <- 'C:/Users/Q/Dropbox/projects/nasabiodiv' # Can be changed to /mnt/research/nasabio/data/bbs to get data from HPCC
fpfig <- 'C:/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_diversity_maps/' # Directory to save images

bd <- read.csv(file.path(fpdata, 'bbs_beta_byroute.csv'), stringsAsFactors = FALSE) # older beta-div, with td pd and fd
bdpart <- read.csv(file.path(fpdata, 'bbs_betapart_byroute.csv'), stringsAsFactors = FALSE) # newer beta-div, with partitioning but no fd
ad <- read.csv(file.path(fpdata, 'bbs_alphadiv.csv'), stringsAsFactors = FALSE) # alpha
gd <- read.csv(file.path(fpdata, 'bbs_gammadiv.csv'), stringsAsFactors = FALSE) # gamma
ed <- read.csv(file.path(fpdata, 'bbs_geodiversity_stats.csv'), stringsAsFactors = FALSE) # all geodiversity

# Calculate means of all routes across years 1997-2016

library(dplyr)

bdpart_yrmeans <- bdpart %>%
  filter(year >= 1997) %>%
  group_by(rteNo, lon, lat, radius, diversity, partition) %>%
  summarize(beta = mean(na.omit(beta)))

# Bivariate plots ---------------------------------------------------------


# Maps --------------------------------------------------------------------

library(cowplot)
source('figs/bbs_map_drawing_fns.r')

radii <- c(50000, 100000, 200000, 300000)

colscalebeta <- scale_colour_gradientn(name = 'Taxonomic\nbeta-diversity', colours = colorRampPalette(colors=RColorBrewer::brewer.pal(9, 'YlOrRd'), bias = 0.2)(9), breaks = c(0,.5,1), limits=c(0,1))


beta_td_total <- bdpart_yrmeans %>%
  filter(diversity == 'taxonomic', partition == 'total', radius %in% radii, !is.na(beta))
  
draw_bbs_map_wrap(dat = beta_td_total, zvar = 'beta', rad = radii, 
                  colscale = colscalebeta,
                  maptitle = 'BBS beta-diversity, taxonomic, incidence-based, total',
                  fp = fpfig, fname = 'beta_div_tax_total.png')
