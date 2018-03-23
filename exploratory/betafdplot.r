# Inspect beta functional diversity values. See where there is a discrepancy

library(dplyr)
library(tidyr)
library(reshape2)

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'
bbsbio <- read.csv(file.path(fp, 'bbs_allbio_wide.csv'), stringsAsFactors = FALSE)

library(cowplot)
p <- ggplot(bbsbio) + theme_bw()

p + geom_histogram(aes(x = beta_fd_pairwise_pa_100)) # Raw
p + geom_histogram(aes(x = beta_fd_pairwise_pa_z_100)) # SES

table(bbsbio$beta_fd_pairwise_pa_z_100 < -1000) # Only one outlier
table(bbsbio$beta_fd_pairwise_pa_z_100 < -100)
p + geom_histogram(aes(x = beta_fd_pairwise_pa_z_100)) + scale_x_continuous(limits = c(-20,10))

p + geom_point(aes(x = beta_fd_pairwise_pa_z_100, y = beta_td_pairwise_pa_100)) # Only a few very low values, and they're all only moderate in taxonomic div.
p + geom_point(aes(x = beta_fd_pairwise_pa_z_100, y = beta_pd_pairwise_pa_z_100)) # Also only moderate with phylogenetic div.
p + geom_point(aes(x = beta_fd_pairwise_pa_z_100, y = alpha_richness_100)) # Same with richness
