# Script to generate qsub bash commands to do all FIA geodiversity calculations.

x <- read.csv('spatial_data_extraction/geodiv_stat_table_separate.csv', stringsAsFactors = FALSE)

library(dplyr)

qsub_string <- function(taxon, var, mem, start, end, walltime='4:00:00') paste0('qsub geoextract.sh -N ', var, '_', taxon, ' -v taxon=', taxon, ',geovar=', var, ' -l walltime=', walltime, ',mem=', mem, 'gb -t ', start, '-', end)

qsub_calls <- x %>%
  rowwise %>%
  do(qsubs = sapply(seq(1, .$N.slices.fiaall, by = 250), function(x) qsub_string(taxon = 'fia', var = .$variable.id, mem = .$RAM, start = x, end = min(x+249, .$N.slices.fiaall))))

write.table(unlist(qsub_calls$qsubs), file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/geo_qsub_all.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

# Do this for the TRI rasters too. (BBS and FIA)

x <- read.csv('spatial_data_extraction/geodiv_stat_table_separate.csv', stringsAsFactors = FALSE)
x <- x[-(1:13),]

library(dplyr)

qsub_string <- function(taxon, var, mem, start, end, walltime='4:00:00') paste0('qsub geoextract.sh -N ', var, '_', taxon, ' -v taxon=', taxon, ',geovar=', var, ' -l walltime=', walltime, ',mem=', mem, 'gb -t ', start, '-', end)

qsub_calls <- x %>%
  rowwise %>%
  do(qsubs = sapply(seq(1, .$N.slices.fiaall, by = 250), function(x) qsub_string(taxon = 'fia', var = .$variable.id, mem = .$RAM, start = x, end = min(x+249, .$N.slices.fiaall))))

write.table(unlist(qsub_calls$qsubs), file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/tri_qsub.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

qsub_calls_bbs <- x %>%
  rowwise %>%
  do(qsubs = sapply(seq(1, .$N.slices.bbs, by = 250), function(x) qsub_string(taxon = 'bbs', var = .$variable.id, mem = .$RAM, start = x, end = min(x+249, .$N.slices.bbs))))

write.table(unlist(qsub_calls_bbs$qsubs), file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/tri_qsub.txt', quote = FALSE, col.names = FALSE, row.names = FALSE, append = TRUE)
