# Use stratified sampling to generate a subset of the BBS data.
# Also ensure that those routes have all data.

r <- 100e3 # minimum distance in m to separate all plots.

source('stats/spatial_fns.r')
source('stats/SRS_iterative.r')

fp <- 'C:/Users/Q/Dropbox/projects/nasabiodiv'

bbsbioall <- read.csv(file.path(fp, 'bbs_allbio_wide.csv'), stringsAsFactors = FALSE)
bbsgeoall <- read.csv(file.path(fp, 'bbs_allgeo_wide.csv'), stringsAsFactors = FALSE)

# Find route names and coordinates that have no missing data
# Nothing is in MNTDfunc so have to get rid
bbsbioall <- bbsbioall[,-grep('MNTDfunc', names(bbsbioall))]

bbsrtes <- bbsbioall[complete.cases(bbsbioall), 1:5]

library(sp)

# Create spatial info and throw out the plots too close to other countries.
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
bbssp <- SpatialPoints(coords = bbsrtes[,c('lon_aea','lat_aea')], proj4string = CRS(aea_crs))

# Throw out the ones that are within 100 km of Mexico or Canada (about 200-300)
bbs_border100 <- flag_coast_plots(focal_points = bbssp, radius = r, border_countries = c('Canada', 'Mexico'))

bbsrtes <- cbind(bbsrtes, bbs_border100)
bbs_noedge <- bbssp[!bbsrtes$is_edge]
bbsrtes_noedge <- subset(bbsrtes, !is_edge)

set.seed(1001)

sample_idx <- SRS_iterative_N1(focal_points = bbs_noedge, radius = r, n = 200) # Returns 153 points.

rteids <- bbsrtes_noedge$rteNo[sample_idx]

# Put the biodiversity and geodiversity data from these routes together into a single data frame.

bbsbiosub <- subset(bbsbioall, rteNo %in% rteids)
bbsgeosub <- subset(bbsgeoall, rteNo %in% rteids)

bbssub <- left_join(bbsbiosub, bbsgeosub)

write.csv(bbssub, file = 'C:/Users/Q/google_drive/NASABiodiversityWG/SampleData/bbs_subset_wideform.csv', row.names = FALSE)


#################
# Pull FIA (run on remote)

r <- 100e3 # minimum distance in m to separate all plots.

source('/mnt/research/nasabio/code/spatial_fns.r')
source('/mnt/research/nasabio/code/SRS_iterative.r')

fp <- '/mnt/research/nasabio/data/fia'

fiabioall <- read.csv(file.path(fp, 'fia_allbio_wide.csv'), stringsAsFactors = FALSE)

# Load geodiv later when needed.

# Find route names and coordinates that have no missing data

fiaplots <- fiabioall[complete.cases(fiabioall), 'PLT_CN']

source('/mnt/research/nasabio/code/loadfiaall.r')

library(sp)

goodfiacoords <- fiacoords[match(fiaplots, fiacoords$PLT_CN),]

# Create spatial info and throw out the plots too close to other countries.
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
fia_aea <- spTransform(SpatialPoints(coords=goodfiacoords[,c('lon','lat')], proj4string = CRS('+proj=longlat')), CRSobj = CRS(aea_crs))

# Throw out the ones that are within 100 km of Mexico or Canada (about 200-300)
fia_border100 <- flag_coast_plots(focal_points = fia_aea, radius = r, border_countries = c('Canada', 'Mexico'))

goodfiacoords <- cbind(goodfiacoords, fia_border100)
fia_noedge <- fia_aea[!goodfiacoords$is_edge]
goodfiacoords_noedge <- subset(goodfiacoords, !is_edge)

set.seed(1002)

sample_idx <- SRS_iterative_N1(focal_points = fia_noedge, radius = r, n = 200)  # Returns 142

pltcns <- goodfiacoords_noedge$PLT_CN[sample_idx]

# Put the biodiversity and geodiversity data from these routes together into a single data frame.

fiabiosub <- subset(fiabioall, PLT_CN %in% pltcns)

# Sequentially load each geo data frame, take a subset, and add to what we have.

geo_list <- list()

fiageopt <- read.csv(file.path(fp, 'fia_geo_by_point.csv'), stringsAsFactors = FALSE)
geo_list[[1]] <- subset(fiageopt, PLT_CN %in% pltcns)
rm(fiageopt)

for (i in dir(file.path(fp, 'geodiv'), pattern = 'wide')) {
  print(i)
  x <- read.csv(file.path(fp, 'geodiv', i))
  geo_list[[length(geo_list) + 1]] <- subset(x, PLT_CN %in% pltcns)
  rm(x)
}

fiasub <- Reduce(left_join, geo_list)
fiasub <- left_join(fiabiosub, fiasub)
write.csv(fiasub, file = '/mnt/research/nasabio/data/fia/fia_subset_wideform.csv', row.names = FALSE)


###
# Get the columns in the right place
fia_dat <- fia_dat %>% select(PLT_CN, HUC4, lon,lat,lon_aea,lat_aea,everything())
bbs_dat <- bbs_dat %>% select(rteNo, HUC4, lon,lat,lon_aea,lat_aea,everything())
write.csv(fia_dat, file = 'C:/Users/Q/google_drive/NASABiodiversityWG/SampleData/fia_subset_wideform.csv', row.names = FALSE)
write.csv(bbs_dat, file = 'C:/Users/Q/google_drive/NASABiodiversityWG/SampleData/bbs_subset_wideform.csv', row.names = FALSE)
