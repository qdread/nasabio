# Convert .RData objects for model fitting into .csv files (better for archiving with manuscript)
# QDR/nasabioxgeo/04 Oct 2018

setwd('/mnt/research/nasabio/temp')

load('bbs_spatial_mm_dat_50k.RData')

dimnames(tnc_bin)[[2]] <- dimnames(tnc_bin)[[1]]

write.csv(tnc_bin, 'tnc_adjacencymatrix.csv')

# Check to make sure the two are the same
tnc2 <- read.csv('tnc_adjacencymatrix.csv', row.names = 1)

tnc2 <- as.matrix(tnc2)
dimnames(tnc2)[[2]] <- NULL

write.csv(bbsbio, 'bbs_biodiversity.csv', row.names = FALSE)
write.csv(bbsgeo, 'bbs_geodiversity.csv', row.names = FALSE)

load('fia_spatial_mm_dat_50k.RData')
write.csv(fiabio, 'fia_biodiversity.csv', row.names = FALSE)
write.csv(fiageo, 'fia_geodiversity.csv', row.names = FALSE)

