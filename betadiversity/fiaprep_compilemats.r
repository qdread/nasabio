# compile matrices for fia

fp <- '/mnt/research/nasabio/data/fia/mats/'

for (r in c(100, 150, 200, 300) * 1000) {
  mats_list <- list()
  for (i in 1:100) {
    load(paste0(fp, 'unfuzzedmat_', as.character(as.integer(r)), '_', i, '.r'))
    mats_list <- c(mats_list, all_mats)
  }
  all_mats <- mats_list
  save(all_mats, file = paste0(fp, 'unfuzzedmat_', as.character(as.integer(r)), '.r'))
}