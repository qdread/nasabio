# Paste together sbatch calls for not completed beta diversity jobs

n <- 119177
h <- 12

notdone <- which(!file.exists(paste0('/mnt/research/nasabio/data/fia/diversity/usa2018/beta_',1:n,'.r')))

#new_calls <- paste0('sbatch --time=12:00:00 --array=', ifelse(notdone %% 1000 == 0, 1000, notdone %% 1000), ' --export=N1000=', ifelse(notdone %% 1000 == 0, floor(notdone/1000) - 1, floor(notdone/1000)), ' fiabeta.sb')

# Better version to make only one sbatch for each group of 1000.
# Better yet break it into groups of 100.
hundreds <- floor(notdone/100) - ifelse(notdone %% 100==0, 1, 0)
new_calls <- sapply(unique(hundreds), function(i) {
	notdone100 <- notdone[hundreds == i] %% 1000
	paste0('sbatch --time=', h, ':00:00 --array=', paste(ifelse(notdone100 == 0, 1000, notdone100), collapse = ','), ' --export=N1000=', floor(notdone[hundreds == i][1] / 1000), ' fiabeta.sb')
})

write.table(new_calls, quote = FALSE, row.names = FALSE)
