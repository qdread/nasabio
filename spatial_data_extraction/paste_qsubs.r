# Script to generate qsub bash commands to do all FIA geodiversity calculations.


# Edit for Jan 18: beta diversity
starts <- seq(1, 135174, by = 250)
isbig <- starts >= 100000

qsub_calls <- paste('qsub fiabd.sh -v isbig=', ifelse(isbig,'yes','no'), ' -t ', ifelse(isbig, as.integer(starts-100000), as.integer(starts)), '-', ifelse(isbig, as.integer(starts-100000+249), as.integer(starts+249)), sep='')
write.table(qsub_calls, file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/bd_qsub.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

# Edited 27 Nov 2018: beta diversity for FIA using SLURM
# Edited 28 Nov 2018: break into groups of 500 instead of 1000 so that the jobs don't get stuck as easily.
# Edited 12 Dec 2018: try 250 instead.

n_plots <- 119177
n_calls <- ceiling(n_plots/250)
sb_calls <- rep(c('sbatch --array=1-250', 'sbatch --array=251-500', 'sbatch --array=501-750', 'sbatch --array=751-1000'), times = floor(n_plots/1000))
sb_calls <- paste(sb_calls, ' --export=N1000=', rep(0:(floor(n_plots/1000) - 1), each = 4), ' fiabeta.sb', sep = '')
sb_calls <- c(sb_calls, 'sbatch --array=1-177 --export=N1000=119 fiabeta.sb')

write.table(sb_calls, file = '~/Dropbox/projects/nasabiodiv/code/bd_qsub.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

# BBS submissions
bbs_sb_calls <- rep(c('sbatch --array=1-250', 'sbatch --array=251-500', 'sbatch --array=501-750', 'sbatch --array=751-1000'), times = 2)
bbs_sb_calls <- paste(bbs_sb_calls, rep(c('--export=N1000=0 bbsbeta.sb', '--export=N1000=1 bbsbeta.sb'), each = 4))
write.table(bbs_sb_calls, quote = FALSE, col.names = FALSE, row.names = FALSE)


# Elevation GDAL-only qsubs
starts <- seq(1, 10000, by = 250)
vars <- c('elevation', 'slope', 'roughness', 'tri')

qsub_calls <- apply(expand.grid(starts, vars), 1, function(x) paste0('qsub elevextract.sh -N ', x[2], ' -v geovar=', x[2], ' -l mem=1gb -t ', x[1], '-', as.character(as.integer(min(as.numeric(x[1])+249, 135174)))))
write.table(qsub_calls, file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/elev_qsub.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)


###########################################
# 11 Jan: GDAL only. All for both BBS and FIA. Delete the ones that are already done.
# Edit 19 Jan: add aspect sin and cos since we have now calculated that!

x <- read.csv('spatial_data_extraction/geodiv_table_for_gdal.csv', stringsAsFactors = FALSE)

library(dplyr)

qsub_string <- function(taxon, var, mem, tmpmem, start, end, walltime='4:00:00') paste0('qsub geoextract.sh -N ', var, '_', taxon, ' -v taxon=', taxon, ',geovar=', var, ' -l walltime=', walltime, ',mem=', mem, 'gb,file=', tmpmem, 'gb -t ', start, '-', end)

qsub_calls_fia <- x %>%
  rowwise %>%
  do(qsubs = sapply(seq(1, .$N.slices.fiaall, by = 250), function(x) qsub_string(taxon = 'fia', var = .$variable.id, mem = 2, tmpmem = 4, start = as.character(as.integer(x)), end = as.character(as.integer(min(x+249, .$N.slices.fiaall))))))

qsub_calls_bbs <- x %>%
  rowwise %>%
  do(qsubs = sapply(seq(1, .$N.slices.bbs, by = 250), function(x) qsub_string(taxon = 'bbs', var = .$variable.id, mem = 2, tmpmem = 4, start = as.character(as.integer(x)), end = as.character(as.integer(min(x+249, .$N.slices.bbs))))))

# Remove whatever has already been completed.
#qsub_calls_fia <- qsub_calls_fia[!x$variable.id %in% c('elevation','slope','roughness','tri'), ]
#qsub_calls_bbs <- qsub_calls_bbs[!x$variable.id %in% c('aspect_sin', 'aspect_cos'), ]

write.table(unlist(qsub_calls_bbs$qsubs), file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/geo_qsub_all_bbs.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

write.table(c(unlist(qsub_calls_fia$qsubs), unlist(qsub_calls_bbs$qsubs)), file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/geo_qsub_all.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)


#############################################
# 28 Jun: Do one job per geo variable, to get 1 km radius stats only (not done in previous, unfortunately)

qsub_string_1job <- function(taxon, var, mem, tmpmem, walltime='4:00:00') paste0('qsub geoextract1km.sh -N ', var, '_', taxon, ' -v taxon=', taxon, ',geovar=', var, ' -l walltime=', walltime, ',mem=', mem, 'gb,file=', tmpmem, 'gb')

qsub_calls_fia <- x %>%
  rowwise %>%
  transmute(qsubs = qsub_string_1job(taxon = 'fia', var = variable.id, mem = 2, tmpmem = 4))

qsub_calls_bbs <- x %>%
  rowwise %>%
  transmute(qsubs = qsub_string_1job(taxon = 'bbs', var = variable.id, mem = 2, tmpmem = 4))
  
write.table(c(qsub_calls_fia$qsubs, qsub_calls_bbs$qsubs), file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/code/geo_qsub_all_1km.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
