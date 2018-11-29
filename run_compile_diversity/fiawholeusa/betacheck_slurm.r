# Paste together sbatch calls for not completed beta diversity jobs

notdone <- which(!file.exists(paste0('/mnt/research/nasabio/data/fia/diversity/usa2018/beta_',1:135032,'.r')))

new_calls <- paste0('sbatch --array=', ifelse(notdone %% 1000 == 0, 1000, notdone %% 1000), ' --export=N1000=', ifelse(notdone %% 1000 == 0, floor(notdone/1000) - 1, floor(notdone/1000)), ' fiabeta.sb')

write.table(new_calls, quote = FALSE, row.names = FALSE)