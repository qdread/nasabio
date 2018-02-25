# Create metadata table for column names

cnames <- union(names(bbs_dat), names(fia_dat))

metadata <- data.frame(name = cnames,
                       variable_type = 'geo',
                       diversity_type = NA,
                       diversity_level = NA,
                       variable_category = NA,
                       variable_subcategory = sapply(strsplit(cnames, '_'), '[', 1),
                       summary_statistic = NA,
                       resolution = NA,
                       radius = 0,
                       stringsAsFactors = FALSE)

metadata$variable_type[grep('alpha|beta|gamma', metadata$name)] <- 'bio'
metadata$variable_type[metadata$name %in% c('rteNo','PLT_CN','lat','lon','lat_aea','lon_aea','HUC4')] <- 'id'
for (s in c('mean','sd','tri','roughness')) metadata$summary_statistic[grep(paste0('_',s), metadata$name)] <- s
metadata$summary_statistic[grep('mode', metadata$name)] <- 'mode'
metadata$summary_statistic[grep('richness_geodiv', metadata$name)] <- 'richness'
metadata$summary_statistic[grep('diversity_geodiv', metadata$name)] <- 'diversity'

for (r in c(5,10,20,30,40,50,75,100,200,300,400,500)) {
  metadata$radius[grep(paste0('_',r,'_'), metadata$name)] <- r
  metadata$radius[grep(paste0('_',r,'$'), metadata$name)] <- r
}
for (d in c('alpha', 'beta', 'gamma')) metadata$diversity_level[grep(d, metadata$name)] <- d
tax_names <- c('richness','shannon','evenness', '_td_')
phy_names <- c('_MPD_', '_MNTD_','_pd_')
func_names <- c('MPDfunc', 'MNTDfunc','_fd_')

for (n in tax_names) metadata$diversity_type[metadata$variable_type == 'bio' & grepl(n, metadata$name)] <- 'taxonomic'
for (n in phy_names) metadata$diversity_type[metadata$variable_type == 'bio' & grepl(n, metadata$name)] <- 'phylogenetic'
for (n in func_names) metadata$diversity_type[metadata$variable_type == 'bio' & grepl(n, metadata$name)] <- 'functional'

metadata$resolution[metadata$variable_type == 'geo'] <- 'native' 
metadata$resolution[metadata$variable_type == 'geo' & grepl('_5k_', metadata$name)] <- 'standard'

metadata$variable_category[metadata$variable_type == 'geo'] <- 'bioclim'
metadata$variable_category[grep('cloud', metadata$name)] <- 'biocloud'
metadata$variable_category[grep('elevation|aspect|slope', metadata$name)] <- 'topography'
metadata$variable_category[grep('dhi', metadata$name)] <- 'DHI'
metadata$variable_category[grep('soil|geo', metadata$name)] <- 'geology'
metadata$variable_category[grep('light|foot', metadata$name)] <- 'humans'

metadata$which_dataset <- ifelse(metadata$name %in% names(bbs_dat) & metadata$name %in% names(fia_dat), 'both',
                                 ifelse(metadata$name %in% names(bbs_dat), 'BBS', 'FIA'))

library(dplyr)
metadata <- metadata %>%
  mutate(variable_type = factor(variable_type, levels = c('id','bio','geo'))) %>%
  arrange(variable_type, diversity_level, diversity_type, variable_category, variable_subcategory, summary_statistic, resolution, radius)
metadata <- metadata[c(6,7,1:5, 8:nrow(metadata)),]
write.csv(metadata, file = 'C:/Users/Q/google_drive/NASABiodiversityWG/SampleData/metadata.csv', row.names = FALSE)
