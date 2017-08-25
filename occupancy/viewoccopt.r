# Read occupancy CSV files and see if they have good estimates of detection probability

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv/occupancy'

occ_output <- list()

for (i in 1:20) {
  occ_output[[i]] <- read.csv(file.path(fp, paste0('occ_ranef_opt', i+1996, '.csv')), comment.char = '#')
}

occ_output <- cbind(year = 1997:2016, do.call('rbind', occ_output))

# I think it completely failed.