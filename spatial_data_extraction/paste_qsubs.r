# Script to generate qsub bash commands to do all FIA geodiversity calculations.

x <- read.csv('spatial_data_extraction/geodiv_stat_table_separate.csv', stringsAsFactors = FALSE)

library(dplyr)

qsub_string <- function(taxon, var, mem, start, end, walltime='4:00:00') paste0('qsub geoextract.sh -N ', var, '_', taxon, ' -v taxon=', taxon, ',geovar=', var, ' -l walltime=', walltime, ',mem=', mem, 'gb -t ', start, '-', end)

qsub_calls <- x %>%
  rowwise %>%
  do(qsubs = sapply(seq(1, .$N.slices.fiaall, by = 250), function(x) qsub_string(taxon = 'fia', var = .$variable.id, mem = .$RAM, start = as.character(as.integer(x)), end = as.character(as.integer(min(x+249, .$N.slices.fiaall))))))

write.table(unlist(qsub_calls$qsubs), file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/geo_qsub_all.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

# Do this for the TRI rasters too. (BBS and FIA)

x <- read.csv('spatial_data_extraction/geodiv_stat_table_separate.csv', stringsAsFactors = FALSE)
x <- x[-(1:13),]

library(dplyr)

qsub_string <- function(taxon, var, mem, start, end, walltime='4:00:00') paste0('qsub geoextract.sh -N ', var, '_', taxon, ' -v taxon=', taxon, ',geovar=', var, ' -l walltime=', walltime, ',mem=', mem, 'gb -t ', start, '-', end)

qsub_calls <- x %>%
  rowwise %>%
  do(qsubs = sapply(seq(1, .$N.slices.fiaall, by = 250), function(x) qsub_string(taxon = 'fia', var = .$variable.id, mem = .$RAM, start = as.character(as.integer(x)), end = as.character(as.integer(min(x+249, .$N.slices.fiaall))))))

write.table(unlist(qsub_calls$qsubs), file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/tri_qsub.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

qsub_calls_bbs <- x %>%
  rowwise %>%
  do(qsubs = sapply(seq(1, .$N.slices.bbs, by = 250), function(x) qsub_string(taxon = 'bbs', var = .$variable.id, mem = .$RAM, start = as.character(as.integer(x)), end = as.character(as.integer(min(x+249, .$N.slices.bbs))))))

write.table(unlist(qsub_calls_bbs$qsubs), file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/tri_qsub.txt', quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)

# Beta-diversity qsubs

starts <- seq(1, 135174, by = 250)
qsub_calls <- paste('qsub fiabd.sh -t ', starts, '-', starts+249, sep='')
write.table(qsub_calls, file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/bd_qsub.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

# Elevation GDAL-only qsubs
starts <- seq(1, 10000, by = 250)
vars <- c('elevation', 'slope', 'roughness', 'tri')

qsub_calls <- apply(expand.grid(starts, vars), 1, function(x) paste0('qsub elevextract.sh -N ', x[2], ' -v geovar=', x[2], ' -l mem=1gb -t ', x[1], '-', as.character(as.integer(min(as.numeric(x[1])+249, 135174)))))
write.table(qsub_calls, file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/elev_qsub.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)


###########################################
# 11 Jan: GDAL only. All for both BBS and FIA. Delete the ones that are already done.

x <- read.csv('spatial_data_extraction/geodiv_table_for_gdal.csv', stringsAsFactors = FALSE)

library(dplyr)

qsub_string <- function(taxon, var, mem, start, end, walltime='4:00:00') paste0('qsub geoextract.sh -N ', var, '_', taxon, ' -v taxon=', taxon, ',geovar=', var, ' -l walltime=', walltime, ',mem=', mem, 'gb -t ', start, '-', end)

qsub_calls_fia <- x %>%
  rowwise %>%
  do(qsubs = sapply(seq(1, .$N.slices.fiaall, by = 250), function(x) qsub_string(taxon = 'fia', var = .$variable.id, mem = 1, start = as.character(as.integer(x)), end = as.character(as.integer(min(x+249, .$N.slices.fiaall))))))

qsub_calls_bbs <- x %>%
  rowwise %>%
  do(qsubs = sapply(seq(1, .$N.slices.fiaall, by = 250), function(x) qsub_string(taxon = 'bbs', var = .$variable.id, mem = 1, start = as.character(as.integer(x)), end = as.character(as.integer(min(x+249, .$N.slices.bbs))))))

qsub_calls_fia <- qsub_calls_fia[, !x$variable.id %in% c('elevation','slope','roughness','tri','')]

write.table(unlist(qsub_calls$qsubs), file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/elev_qsub_all.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)