# Black and white figures for ESA 2019 powerpoint
# QDR/NASABioxgeo/5 June 2019

# Modified 19 June 2019: make the points and lines bigger and/or thicker.

# Load and combine data ---------------------------------------------------

fp <- '~/Dropbox/projects/nasabiodiv/modelfits' # Local

model_coef <- read.csv(file.path(fp, 'multivariate_spatial_coef.csv'), stringsAsFactors = FALSE)
model_pred <- read.csv(file.path(fp, 'multivariate_spatial_pred.csv'), stringsAsFactors = FALSE)
model_rmse <- read.csv(file.path(fp, 'multivariate_spatial_rmse.csv'), stringsAsFactors = FALSE)
model_r2 <- read.csv(file.path(fp, 'multivariate_spatial_r2.csv'), stringsAsFactors = FALSE)
model_coef_var <- read.csv(file.path(fp, 'multivariate_spatial_coef_variation_corrected.csv'), stringsAsFactors = FALSE)
model_waic <- read.csv(file.path(fp, 'multivariate_spatial_waic.csv'), stringsAsFactors = FALSE)
kfold_rmse <- read.csv(file.path(fp, 'multivariate_kfold_rmse.csv'), stringsAsFactors = FALSE)


library(dplyr)
library(ggplot2)
library(reshape2)
library(purrr)

prednames50 <- c('elevation_5k_tri_50_mean', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'dhi_gpp_5k_tri_50_mean')
geo_names <- c('elevation diversity','temperature mean','geol. age diversity','soil diversity','precip. mean','GPP diversity')
geo_names_order <- c('temperature mean', 'precip. mean', 'elevation diversity', 'GPP diversity', 'geol. age diversity', 'soil diversity')

bio_titles <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'beta PD', 'gamma PD', 'alpha FD', 'beta FD', 'gamma FD')
bio_names <- c("alpha_richness", "beta_td_sorensen_pa", "gamma_richness",
               "alpha_phy_pa", "beta_phy_pa", "gamma_phy_pa", 
               "alpha_func_pa", "beta_func_pa", "gamma_func_pa")

twocolors <- c('black', 'gray60')

# Combine full-model and k-fold RMSEs.
# Include RMSE from each fold so we can see variability due to folds.
# Edited 02 May 2019: correct the data wrangling code because the output of the new k-fold is slightly different.

all_rmse <- kfold_rmse %>%
  rename_if(is.numeric, ~ paste0('kfold_', .x)) %>%
  right_join(model_rmse) %>%
  left_join(model_r2 %>% select(-fold) %>% rename(r2 = Estimate, r2_error = Est.Error, r2_q025 = Q2.5, r2_q975 = Q97.5)) %>%
  mutate(response = factor(bio_titles[match(response, bio_names)], levels = bio_titles),
         flavor = map_chr(strsplit(as.character(response), ' '), 2) %>%
           factor(levels = c('TD','PD','FD'), labels = c('taxonomic', 'phylogenetic', 'functional')))
  

# Reshape coefficient plot and relabel it
all_coef <- model_coef %>%
  filter(effect == 'fixed', !parameter %in% 'Intercept') %>%
  dcast(taxon + rv + model + response + parameter ~ stat) %>%
  mutate(predictor = factor(geo_names[match(parameter, prednames50)], levels = geo_names_order),
         response = factor(bio_titles[match(response, bio_names)], levels = bio_titles),
         flavor = map_chr(strsplit(as.character(response), ' '), 2) %>%
           factor(levels = c('TD','PD','FD'), labels = c('taxonomic', 'phylogenetic', 'functional')))

# Relabel data frame of spatial variability metrics  
model_coef_var <- model_coef_var %>%
  mutate(predictor = factor(geo_names[match(parameter, prednames50)], levels = geo_names_order),
         response = factor(bio_titles[match(response, gsub('_', '', bio_names))], levels = bio_titles),
         flavor = map_chr(strsplit(as.character(response), ' '), 2) %>%
           factor(levels = c('TD','PD','FD'), labels = c('taxonomic', 'phylogenetic', 'functional'))) %>%
  rename(coef_var = Estimate)

# Coefficient plots -------------------------------------------------

fpfig <- '~/google_drive/NASABiodiversityWG/Conferences/ESA2019/talkimgs' 
source('~/Documents/R/theme_black.R')
# Add some color to indicate which ones' credible intervals are not zero
# Also shade the climate mean region with a gray rectangle
# Only include taxonomic diversity, and make it black and white

coefdat_bbs <- all_coef %>%
  filter(taxon == 'bbs', flavor == 'taxonomic') %>%
  mutate(nonzero = Q2.5 > 0 | Q97.5 < 0,
         rv = paste0(rv, '-diversity'))
coefplot_bbs <- ggplot(coefdat_bbs %>% filter(model=='full')) +
  geom_rect(xmin=0, xmax=2.5, ymin=-Inf, ymax=Inf, fill = 'gray50', alpha = 0.05) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(x = predictor, ymin = Q2.5, ymax = Q97.5, color = nonzero), width = 0, size = 1) +
  geom_point(aes(x = predictor, y = Estimate, color = nonzero, size = nonzero)) +
  facet_grid(. ~ rv, labeller = label_parsed) +
  scale_y_continuous(name = 'coefficient estimate', limits = c(-0.73, 0.73), expand = c(0,0)) +
  scale_color_manual(values = c('white', 'red')) +
  scale_size_manual(values = c(2, 3)) +
  theme_black() +
  theme(strip.background = element_rect(fill=NA, color = 'white'),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')

coefdat_fia <- all_coef %>%
  filter(taxon == 'fia', flavor == 'taxonomic') %>%
  mutate(nonzero = Q2.5 > 0 | Q97.5 < 0,
         rv = paste0(rv, '-diversity'))
coefplot_fia <- ggplot(coefdat_fia %>% filter(model=='full')) +
  geom_rect(xmin=0, xmax=2.5, ymin=-Inf, ymax=Inf, fill = 'gray50', alpha = 0.05) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'slateblue', size = 1) +
  geom_errorbar(aes(x = predictor, ymin = Q2.5, ymax = Q97.5, color = nonzero), width = 0, size = 1) +
  geom_point(aes(x = predictor, y = Estimate, color = nonzero, size = nonzero)) +
  facet_grid(. ~ rv, labeller = label_parsed) +
  scale_y_continuous(name = 'coefficient estimate', limits = c(-0.73, 0.73), expand = c(0,0)) +
  scale_color_manual(values = c('white', 'red')) +
  scale_size_manual(values = c(2, 3)) +
  theme_black() +
  theme(strip.background = element_rect(fill=NA, color = 'white'),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'none')

# Edit 6 June: sideways coefficient plot.
coef_bbs_sideways <- coefplot_bbs + 
  coord_flip() +
  theme(axis.text.x = element_text(angle=0, hjust=0.5)) 
coef_fia_sideways <- coefplot_fia + 
  coord_flip() +
  theme(axis.text.x = element_text(angle=0, hjust=0.5)) 

ggsave(file.path(fpfig, 'BBS_coef.png'), coef_bbs_sideways, height = 5*.8, width = 10*.8, dpi = 300)
ggsave(file.path(fpfig, 'FIA_coef.png'), coef_fia_sideways, height = 5*.8, width = 10*.8, dpi = 300)


# Plot showing RMSEs --------------------------------------------------------

pn1 <- position_nudge(x = -0.06, y = 0)
pn2 <- position_nudge(x = 0.06, y = 0)


# Edit 18 June: plot comparing RMSEs and R-squared for the 3 model types
# Edit 08 Aug: add geodiversity-only to this

all_rmse <- all_rmse %>%
  mutate(model = factor(model, levels=c('space','climate','geo','full'), labels=c('space only', 'climate','geodiversity','climate +\ngeodiversity')))

pd <- position_dodge(width = 0.08)



kfold_rmseplot_bymodel_bird <- all_rmse %>% 
  filter(taxon %in% 'bbs', !model %in% 'space only') %>%
  ggplot(aes(x = rv, color = model, group = model)) +
  facet_grid(. ~ flavor, switch = 'x') +
  geom_errorbar(aes(ymin = kfold_RMSE_q025_relative, ymax = kfold_RMSE_q975_relative), width = 0, position = pd) +
  geom_point(aes(y = kfold_RMSE_mean_relative), position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.7), expand = c(0,0), name = 'CV relative root mean squared error') +
  scale_x_discrete(name = 'response') +
  ggtitle('birds') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(0, 'lines'),
        legend.position = 'none')
kfold_rmseplot_bymodel_tree <- all_rmse %>% 
  filter(taxon %in% 'fia', !model %in% 'space only') %>%
  ggplot(aes(x = rv, color = model, group = model)) +
  facet_grid(. ~ flavor, switch = 'x') +
  geom_errorbar(aes(ymin = kfold_RMSE_q025_relative, ymax = kfold_RMSE_q975_relative), width = 0, position = pd) +
  geom_point(aes(y = kfold_RMSE_mean_relative), position = pd) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.7), expand = c(0,0), name = 'CV relative root mean squared error') +
  scale_x_discrete(name = 'response') +
  ggtitle('trees') +
  theme(strip.background = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(0, 'lines'),
        legend.position = c(0.5, 0.2),
        legend.background = element_rect(color = 'black'),
        legend.text = element_text(size = 6.5))

# RMSE plot as 2-way facet

kfold_rmseplot_bymodel_2wayfacet <- all_rmse %>% 
  filter(!model %in% 'space only') %>%
  ggplot(aes(x = rv, color = model, group = model)) +
  facet_grid(taxon ~ flavor, switch = 'x', labeller = labeller(taxon = c(bbs = 'birds', fia = 'trees'))) +
  geom_errorbar(aes(ymin = kfold_RMSE_q025_relative, ymax = kfold_RMSE_q975_relative), width = 0, position = pd, size = 1) +
  geom_point(aes(y = kfold_RMSE_mean_relative), position = pd, size = 2.5) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.7), name = 'CV relative root mean squared error', expand = c(0,0)) +
  scale_x_discrete(name = 'response') +
  theme_black() +
  theme(strip.background = element_blank(),
        strip.placement = 'outside',
        panel.spacing = unit(0, 'lines'),
        legend.position = 'right',
        legend.background = element_rect(color = 'black'),
        legend.direction = 'vertical',
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.text = element_text(size = 8))

ggsave(file.path(fpfig, 'both_lolormse_allmodels_2wayfacet.png'), kfold_rmseplot_bymodel_2wayfacet, height = 6*.8, width = 9*.8, dpi = 400)


# Map of tnc regions ------------------------------------------------------

# Load TNC boundaries
library(sp)
library(rgdal)
library(ggplot2)
library(dplyr)
library(purrr)
library(rgeos)
library(maptools)
fpfig <- '~/google_drive/NASABiodiversityWG/Conferences/ESA2019/talkimgs' 
fpregion <- '~/Dropbox/projects/nasabiodiv/regions'
tnc <- readOGR(dsn = fpregion, layer = 'tnc_terr_ecoregions')
load('~/Documents/R/states_albers.RData')
states <- read.csv('~/Documents/R/states_albers.csv', stringsAsFactors = FALSE)
bbscoords <- read.csv('~/Dropbox/projects/nasabiodiv/bbs_centroids.csv')

# Load coefficients so we can see which regions to keep
fpcoef <- '~/Dropbox/projects/nasabiodiv/modelfits' 
mvcoef <- read.csv(file.path(fpcoef, 'multivariate_spatial_coef.csv'))
regionlist <- unique(mvcoef$region)

aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
tnc <- spTransform(tnc, CRSobj = CRS(aea_crs))

tnc@data <- tnc@data %>%
  mutate(id = rownames(tnc@data), region = as.character(ECODE_NAME))

# Subset out the regions that are outside the US.
tnc <- subset(tnc, region %in% regionlist)

# Clip TNC to US boundaries
goodusabounds <- gUnaryUnion(states_albers)
tncdat <- tnc@data
tnc_clip <- gIntersection(tnc, goodusabounds, byid = TRUE, id = tnc$id)

region_fort <- tnc_clip %>%
  fortify(region = 'id') %>% 
  left_join(tncdat, by = 'id')

blktheme <- theme_bw() + 
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'black'), 
        panel.border = element_blank(), 
        plot.background = element_rect(fill = 'black'), 
        legend.position = c(0.13,0.1), 
        legend.direction = 'horizontal', 
        legend.title = element_blank())

# This is just a map of the region boundaries and the state boundaries.
cols <- RColorBrewer::brewer.pal(5, 'Set2')
testscale <- scale_fill_manual(values = sample(cols, 63, replace = TRUE))

p <- ggplot(region_fort) +
  geom_polygon(aes(x=long, y=lat, group=group, fill=factor(id))) +
  geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.75) +
  geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'gray50', size = 0.75) +
  coord_equal() +
  blktheme

p + testscale

#library(MapColoring)
#tnc_col <- getColoring(tnc_clip)

#hardcoded
tnc_col <- c(1, 2, 4, 3, 2, 1, 1, 4, 3, 2, 2, 2, 3, 2, 1, 1, 2, 3, 3, 2, 
  3, 1, 1, 2, 1, 4, 3, 4, 2, 1, 1, 4, 4, 3, 2, 3, 2, 2, 1, 2, 4, 
  4, 1, 3, 2, 1, 4, 1, 1, 2, 3, 2, 3, 1, 3, 4, 3, 1, 4, 1, 2, 5, 
  4)

cbind(tncdat[,c('id','region')], tnc_col)
# Manually correct so that it needs only 4 colors.
tnc_col_mod <- tnc_col
tnc_col_mod[match(c(745, 770, 732), tncdat$id)] <- c(4,3,3)

# We need to change the drawing order so that the Black Hills shows up on top. 
# Find the ID of the Black Hills and set it to be drawn last.

testscale2 <- scale_fill_manual(values = cols[tnc_col_mod])

p2 <- ggplot(region_fort) +
  geom_polygon(aes(x=long, y=lat, group=group, fill=factor(id))) +
  geom_polygon(data = region_fort %>% filter(grepl('Black Hills', region)), aes(x=long, y=lat, group=group, fill=factor(id))) +
  geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.75) +
  geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'gray50', size = 1) +
  coord_equal() +
  blktheme + theme(legend.position = 'none')

p2 + testscale2
ggsave(file.path(fpfig, 'tncmap.png'), p2 + testscale2, height = 6, width = 9, dpi = 400)


# Map of TNC regions with numbers of sites --------------------------------

# Load TNC boundaries
library(sp)
library(rgdal)
library(ggplot2)
library(dplyr)
library(purrr)
library(rgeos)
library(maptools)
fpfig <- '~/google_drive/NASABiodiversityWG/Conferences/ESA2019/talkimgs' 
fpregion <- '~/Dropbox/projects/nasabiodiv/regions'
tnc <- readOGR(dsn = fpregion, layer = 'tnc_terr_ecoregions')
load('~/Documents/R/states_albers.RData')
states <- read.csv('~/Documents/R/states_albers.csv', stringsAsFactors = FALSE)

# Load bbs and fia ecoregions
bbseco <- read.csv('~/Dropbox/projects/nasabiodiv/bbs_ecoregions.csv', stringsAsFactors = FALSE)
fiaeco <- read.csv('~/Dropbox/projects/nasabiodiv/fia_ecoregions.csv', stringsAsFactors = FALSE)

bbstnc <- bbseco %>%
  filter(!is.na(TNC)) %>%
  group_by(TNC) %>%
  summarize(n_bbs = n())

fiatnc <- fiaeco %>%
  filter(!is.na(TNC)) %>%
  group_by(TNC) %>%
  summarize(n_fia = n())

alltnc <- left_join(bbstnc, fiatnc)

aea_crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
tnc <- spTransform(tnc, CRSobj = CRS(aea_crs))

tnc@data <- tnc@data %>%
  mutate(id = rownames(tnc@data), region = as.character(ECODE_NAME))

# Subset out the regions that are outside the US.
tnc <- subset(tnc, region %in% alltnc$TNC)

# Clip TNC to US boundaries
goodusabounds <- gUnaryUnion(states_albers)
tncdat <- tnc@data
tnc_clip <- gIntersection(tnc, goodusabounds, byid = TRUE, id = tnc$id)

region_fort <- tnc_clip %>%
  fortify(region = 'id') %>% 
  left_join(tncdat, by = 'id') %>%
  left_join(alltnc, by = c('region' = 'TNC'))

blktheme <- theme_bw() + 
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'black'), 
        panel.border = element_blank(), 
        plot.background = element_rect(fill = 'black'), 
        legend.position = c(0.13,0.1), 
        legend.direction = 'horizontal', 
        legend.text = element_text(color = 'white', size = 18),
        legend.background = element_rect(fill = 'transparent'),
        legend.key.width = unit(0.4, 'inches'),
        legend.title = element_blank())

pbbs <- ggplot(region_fort) +
  geom_polygon(aes(x=long, y=lat, group=group, fill=n_bbs)) +
  geom_polygon(data = region_fort %>% filter(grepl('Black Hills', region)), aes(x=long, y=lat, group=group, fill=n_bbs)) +
  geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.75) +
  geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'gray50', size = 1) +
  coord_equal() +
  scale_fill_viridis_c() +
  blktheme 

pfia <- ggplot(region_fort) +
  geom_polygon(aes(x=long, y=lat, group=group, fill=n_fia)) +
  geom_polygon(data = region_fort %>% filter(grepl('Black Hills', region)), aes(x=long, y=lat, group=group, fill=n_fia)) +
  geom_path(aes(x=long, y=lat, group=group), color = 'white', size = 0.75) +
  geom_path(data = states %>% filter(!region %in% c('alaska','hawaii')), aes(x = long, y = lat, group = group), color = 'gray50', size = 1) +
  coord_equal() +
  scale_fill_viridis_c() +
  blktheme 

ggsave(file.path(fpfig, 'tncmap_nbbs.png'), pbbs, height = 6, width = 9, dpi = 400)
ggsave(file.path(fpfig, 'tncmap_nfia.png'), pfia + theme(legend.key.width = unit(0.45, 'inches')) + scale_fill_viridis_c(breaks = c(2500,5000,7500)), height = 6, width = 9, dpi = 400)
