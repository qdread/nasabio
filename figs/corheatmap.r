# Tables and figures of correlations within predictors and within response variables
# QDR/NASABioxgeo/30 May 2018

# Edit 18 Jun: Use new predictor sets

load('C:/Users/Q/Dropbox/projects/nasabiodiv/modelfits/bbs_spatial_mm_dat_50k.RData')
load('C:/Users/Q/Dropbox/projects/nasabiodiv/modelfits/fia_spatial_mm_dat_50k.RData')

prednames50 <- c('elevation_5k_tri_50_mean', 'bio1_5k_50_mean', 'geological_age_5k_50_diversity', 'soil_type_5k_50_diversity', 'bio12_5k_50_mean', 'dhi_gpp_5k_tri_50_mean')
geo_names <- c('elevation diversity','temperature mean','geol. age diversity','soil diversity','precip. mean','GPP diversity')
geo_names_order <- c('temperature mean', 'precip. mean', 'elevation diversity', 'GPP diversity', 'geol. age diversity', 'soil diversity')

bio_titles <- c('alpha TD', 'beta TD', 'gamma TD', 'alpha PD', 'beta PD', 'gamma PD', 'alpha FD', 'beta FD', 'gamma FD')
bio_names <- c("alpha_richness", "beta_td_sorensen_pa", "gamma_richness",
               "alpha_phy_pa", "beta_phy_pa", "gamma_phy_pa", 
               "alpha_func_pa", "beta_func_pa", "gamma_func_pa")

library(dplyr)

fiabio <- fiabio %>%
  select(PLT_CN, !!!bio_names) %>%
  mutate(beta_td_sorensen_pa = qlogis(beta_td_sorensen_pa)) %>%
  setNames(c('PLT_CN', rev(bio_titles)))

bbsbio <- bbsbio %>%
  select(rteNo, !!!bio_names) %>%
  mutate(beta_td_sorensen_pa = qlogis(beta_td_sorensen_pa)) %>%
  setNames(c('PLT_CN', rev(bio_titles)))

fiabio_cor <- round(cor(fiabio[,-1], use = 'pairwise.complete.obs'), 2)
bbsbio_cor <- round(cor(bbsbio[,-1], use = 'pairwise.complete.obs'), 2)

fiageo <- fiageo %>%
  select(PLT_CN, !!!prednames50) 
bbsgeo <- bbsgeo %>%
  select(rteNo, !!!prednames50) 
allgeo <- rbind(fiageo[,-1],bbsgeo[,-1]) %>%
  setNames(geo_names) %>%
  select(!!!rev(geo_names_order))

fiageo_cor <- round(cor(fiageo[,-1], use = 'pairwise.complete.obs'), 2)
bbsgeo_cor <- round(cor(bbsgeo[,-1], use = 'pairwise.complete.obs'), 2)
allgeo_cor <- round(cor(allgeo, use = 'pairwise.complete.obs'), 2)

# Refer here for code on ggplot heatmaps
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

cormat_heatmap <- function(mat) {
  require(reshape2)
  require(ggplot2)
  
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat, diag = TRUE)]<- NA
    return(cormat)
  }
  
  #mat <- reorder_cormat(mat)
  mat <- get_upper_tri(mat)
  
  mat_melt <- melt(mat, na.rm = TRUE)

  ggheatmap <- ggplot(mat_melt, aes(Var1, Var2, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 9, hjust = 1))+
    coord_fixed()
  
  ggheatmap + 
    geom_text(aes(Var1, Var2, label = value), color = "black", size = 2) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(0, 0),
      legend.position = c(0.6, 0.1),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 5, barheight = 0.8,
                                 title.position = "top", title.hjust = 0.5))
  
}

hmfia <- cormat_heatmap(fiabio_cor) + ggtitle('tree biodiversity correlations')
hmbbs <- cormat_heatmap(bbsbio_cor) + theme(legend.position = 'none') + ggtitle('bird biodiversity correlations')
hmgeo <- cormat_heatmap(allgeo_cor) + theme(legend.position = 'none') + ggtitle('geodiversity correlations')

library(gridExtra)
png('C:/Users/Q/google_drive/NASABiodiversityWG/Figures/multivariate_maps_figs/corr_heatmaps.png', height = 8, width = 5, res = 400, units = 'in')
  grid.arrange(hmgeo, hmbbs, hmfia, ncol = 1)
dev.off()
