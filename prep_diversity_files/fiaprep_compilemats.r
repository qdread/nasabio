# compile matrices for fia
# edited 13 Dec for entire USA

fp <- '/mnt/research/nasabio/data/fia/mats/'
radii <- c(5, 10, 20, 50, 75, 100, 150, 200, 300) * 1000
n_tasks <- c(25, 25, 25, 25, 25, 500, 500, 500, 500)

for (j in 1:length(radii)) {
  print(radii[j])
  mats_list <- list()
  for (i in 1:(n_tasks[j])) {
    load(paste0(fp, 'usamat_', as.character(as.integer(radii[j])), '_', i, '.r'))
    mats_list <- c(mats_list, all_mats)
  }
  all_mats <- mats_list
  save(all_mats, file = paste0(fp, 'usamat_', as.character(as.integer(radii[j])), '.r'))
}