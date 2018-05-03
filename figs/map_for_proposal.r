# BBS TNC map formatted for proposal, omitting phylogenetic diversity.
# Requires everything to be loaded from spam_maps.r

maps_bbs_tnc <- coef_bbs %>%
  filter(effect == 'random', ecoregion == 'TNC', parameter != 'Intercept') %>%
  dplyr::select(rv, parameter, region, Estimate) %>%
  rename(ECODE_NAME = region) %>%
  group_by(rv, parameter) %>%
  do(maps = model_map(., rbfill, tnc, states))
elev_maps <- maps_bbs_tnc %>%
  filter(!grepl('phy', rv), parameter == 'elevation_5k_100_sd') %>%
  mutate(bio_title = bio_titles[match(rv, bio_names)])
 
elev_maps <- elev_maps[na.omit(match(bio_titles, elev_maps$bio_title)), ]

whitetheme <-  theme_bw() + 
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(), 
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'gray90', fill = 'gray90'), 
        panel.border = element_blank(), 
        plot.background = element_rect(fill = NA), 
        legend.position = 'none',
        # legend.position = c(0.13,0.1), 
        # legend.direction = 'horizontal', 
        # legend.title = element_blank(),
        # legend.background = element_rect(fill = 'gray90'),
        # legend.text = element_text(color = 'black'),
        plot.title = element_text(color = 'black'))

png('C:/Users/Q/google_drive/NASABiodiversityWG/Figures/bbs_tnc_elevation_slope_map.png', height = 4, width = 8, res = 400, units = 'in')
  grid.arrange(grobs = map2(elev_maps$maps, elev_maps$bio_title, function(p, name) ggplotGrob(p + ggtitle(name) + whitetheme)), nrow = 2)
dev.off()


# Just legend

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

legtheme <- theme(
  legend.position = 'bottom',
  legend.direction = 'horizontal',
  legend.title = element_blank(),
  legend.background = element_rect(fill = NA),
  legend.text = element_blank()
)



legend <- g_legend(elev_maps$maps[[1]] + legtheme) 

library(grid)

png('C:/Users/Q/google_drive/NASABiodiversityWG/Figures/legend.png', height = 1, width = 2, res = 400, units = 'in')
grid.draw(legend) 
dev.off()
