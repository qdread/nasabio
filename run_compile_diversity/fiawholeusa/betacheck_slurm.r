# Paste together sbatch calls for not completed beta diversity jobs

notdone <- which(!file.exists(paste0('/mnt/research/nasabio/data/fia/diversity/usa2018/td100beta_',1:75000,'.r')))

new_calls <- paste0('sbatch --array=', notdone %% 1000, ' --export=N1000=', floor(notdone / 1000), ' fiabeta.sb')

write.table(new_calls, quote = FALSE, row.names = FALSE)