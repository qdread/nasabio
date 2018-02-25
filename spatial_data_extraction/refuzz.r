# Refuzz FIA coordinates so they aren't "too" fuzzy then add them to the generated subset.

fiacoords <- read.csv('~/FIA/FIA10nov/allfia.csv', stringsAsFactors = FALSE)


library(dplyr)
set.seed(42)
fiacoords <- fiacoords %>%
  mutate(latfuzz = ACTUAL_LAT + runif(nrow(fiacoords), -0.005, 0.005),
         lonfuzz = ACTUAL_LON + runif(nrow(fiacoords), -0.005, 0.005))

# Transform the fuzzed coordinates to Albers equal area as well.
library(sp)
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
fia_aea <- spTransform(SpatialPoints(coords=fiacoords[,c('lonfuzz','latfuzz')], proj4string = CRS('+proj=longlat')), CRSobj = CRS(aea_crs))

fiacoords <- mutate(fiacoords, lonfuzz_aea = fia_aea@coords[,1], latfuzz_aea = fia_aea@coords[,2]) %>%
  rename(PLT_CN = CN) %>%
  select(-ACTUAL_LAT, -ACTUAL_LON)

write.csv(fiacoords, '~/FIA/fiafuzzedbyQ.csv', row.names = FALSE)

fia_dat <- read.csv('C:/Users/Q/google_drive/NASABiodiversityWG/SampleData/fia_subset_wideform.csv', stringsAsFactors = FALSE)
fia_dat <- fia_dat %>% 
  left_join(fiacoords) %>%
  setNames(gsub('fuzz','',names(.))) %>%
  select(PLT_CN, lon, lat, lon_aea, lat_aea, HUC4, everything())
write.csv(fia_dat, 'C:/Users/Q/google_drive/NASABiodiversityWG/SampleData/fia_subset_wideform.csv', row.names = FALSE)
