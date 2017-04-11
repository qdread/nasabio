# Compile the beta diversity from fia

fp <- '/mnt/research/nasabio'
load(file.path(fp, 'data/fia/fia_taxonomicbetadiv3.RData')
fia_betadiv3k <- data.frame(beta_td_shannon = fia_shannonbetadiv, beta_td_pairwise = fia_meanpairwisedissim, beta_pd_nearesttaxon_z = fia_phypairwise_z, beta_pd_pairwise_z = fia_phynt_z)
save(fia_betadiv3k, file = file.path(fp, 'data/fia/fia_betadiv3k.r'))

#########################

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/'
load(file.path(fp,'fia_diversitymetrics.RData'))
load(file.path(fp,'fia_betadiv3k.r'))

fia_betadiv3k <- cbind(fiaalbers@data, fia_betadiv3k)
plot_diversity <- left_join(plot_diversity, fia_betadiv3k)
names(plot_diversity)[13:16] <- paste0(names(plot_diversity)[13:16], '3000')