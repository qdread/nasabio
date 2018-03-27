# Generate a single matrix 135174 x 135174 for each of the metrics.

column <- as.numeric(Sys.getenv('PBS_ARRAYID'))

n_plot <- 135174

metric_list <- list()
fp <- '/mnt/research/nasabio/data/fia/diversity/usa'

for (i in 1:n_plot) {
	load(file.path(fp, paste0('beta_', i, '.r')))
	metric_list[[i]] <- beta_div[,column]
	if(i%%1000==0) print(i)
}

library(dplyr)
metric_list <- do.call(rbind, metric_list)
save(metric_list, file = file.path(fp, paste0('matrix_beta_', column, '.r'))
