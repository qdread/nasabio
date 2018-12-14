# Paste together sbatch calls for not completed beta diversity jobs

n <- 119177

notdone <- which(!file.exists(paste0('/mnt/research/nasabio/data/fia/diversity/usa2018/beta_',1:n,'.r')))

#new_calls <- paste0('sbatch --time=12:00:00 --array=', ifelse(notdone %% 1000 == 0, 1000, notdone %% 1000), ' --export=N1000=', ifelse(notdone %% 1000 == 0, floor(notdone/1000) - 1, floor(notdone/1000)), ' fiabeta.sb')

# Better version to make only one sbatch for each group of 1000.
thousands <- floor(notdone/1000) - ifelse(notdone%%1000==0, 1, 0)
new_calls <- sapply(unique(thousands), function(i) {
	notdone1k <- notdone[thousands == i] %% 1000
	paste0('sbatch --time=12:00:00 --array=', paste(ifelse(notdone1k == 0, 1000, notdone1k), collapse = ','), ' --export=N1000=', i, ' fiabeta.sb')
})

write.table(new_calls, quote = FALSE, row.names = FALSE)
