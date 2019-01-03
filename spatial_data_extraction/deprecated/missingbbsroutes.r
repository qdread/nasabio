# Where are the mismatches of BBS routes?

# bbsdf and bbscoords sd be loaded

bbsmeta <- read.csv('~/bbsrtesmetadata.csv',stringsAsFactors=F)

rteNo_birdcounts <- unique(bbsdf$rteNo)
rteNo_shapefile <- unique(bbsmeta$rteno)

table(rteNo_birdcounts %in% rteNo_shapefile)
table(rteNo_shapefile %in% rteNo_birdcounts)

bbscoords$insurvey <- rteNo_shapefile %in% rteNo_birdcounts

bbsdf$hascoords <- bbsdf$rteNo %in% rteNo_birdcounts[rteNo_birdcounts %in% rteNo_shapefile]
state_table <- with(bbsdf, table(statenum, hascoords))

lower48

library(maps)
map('state')
with(bbscoords, points(x = lon, y = lat, col = insurvey + 1))

# get state codes 
statecodes <- read.csv('~/statecodes.csv',stringsAsFactors=FALSE)
state_table[dimnames(state_table)[[1]] %in% statecodes$RegionCode[statecodes$countrynum==840],]

# Find the bbs routes that are in the lower 48 (incl. DC) and that aren't in the shape file.
lower48 <- statecodes$RegionCode[statecodes$countrynum==840 & !statecodes$State %in% c('ALASKA','PUERTO RICO')]
bbsdf$lower48 <- bbsdf$statenum %in% lower48

bbsdf_missing <- bbsdf[!bbsdf$hascoords & bbsdf$lower48, c('statenum', 'Route', 'rteNo')]
bbsdf_missing <- unique(bbsdf_missing)
write.csv(bbsdf_missing, '~/bbsdf_missing.csv', row.names = FALSE)