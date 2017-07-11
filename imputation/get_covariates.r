# Get environmental predictors for the different species locations 

library(raster)
try_spp_loc <- read.csv('/mnt/research/nasabio/data/fia/tree_locations.csv', stringsAsFactors = FALSE)

extract_coords <- with(try_spp_loc, cbind(lon,lat))

srtm_dem <- raster('/mnt/research/ersamlab/shared_data/DEM_SRTM_30m/VRTs/conus_30m_dem.vrt')
try_elevs <- extract(srtm_dem, extract_coords, method = 'simple')

bioclim_stack <- stack('/mnt/research/nasabio/data/bioclim/Bioclim1k/rasterstack/bioclim1k_20170613.vrt')
try_bioclim <- extract(bioclim_stack, extract_coords, method = 'simple')

# NPP


nppnames <- dir('/mnt/research/aquaxterra/DATA/raw_data/MODIS/MOD17A3_gpp_npp', pattern = '*.tif$')
nppextract <- list()

for (i in nppnames) {
  npp_i <- raster(file.path('/mnt/research/aquaxterra/DATA/raw_data/MODIS/MOD17A3_gpp_npp', i))
  nppextract[[length(nppextract) + 1]] <- extract(npp_i, extract_coords, method = 'simple')
  print(i)
}

nppextract <- as.data.frame(do.call('cbind', nppextract))
names(nppextract) <- substr(nppnames, 1, nchar(nppnames)-4)

# Replace 65535 (missing value code) with NA.
nppextract[nppextract == 65535] <- NA

# LAI (convert to AEA for this purpose)
aea_crs <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0'
wgs_crs <- '+proj=longlat +ellps=WGS84 +no_defs'
extract_coords_spatial <- SpatialPoints(extract_coords, proj4string = CRS(wgs_crs))
extract_coords_albers <- spTransform(extract_coords_spatial, CRSobj = CRS(aea_crs))

lainames <- dir('/mnt/research/aquaxterra/DATA/reprojected_data/MODIS/MOD15A2_LAI_FPAR', pattern = 'lai_max.*tif$')
laiextract <- list()

for (i in lainames) {
  lai_i <- raster(file.path('/mnt/research/aquaxterra/DATA/reprojected_data/MODIS/MOD15A2_LAI_FPAR', i))
  laiextract[[length(laiextract) + 1]] <- extract(lai_i, extract_coords_albers, method = 'simple')
  print(i)
}

laiextract <- as.data.frame(do.call('cbind', laiextract))
names(laiextract) <- paste('lai_max', 2001:2011, sep = '_')

# Take median of LAI and NPP over time
lai_median <- apply(laiextract, 1, median, na.rm=T)
npp_median <- apply(nppextract, 1, median, na.rm=T) * 0.1

# Combine and save all covariates

dimnames(try_bioclim)[[2]] <- paste0('bio',1:19)

all_covariates <- cbind(try_spp_loc, lat_aea = extract_coords_albers@coords[,2], lon_aea = extract_coords_albers@coords[,1], elevation = try_elevs, try_bioclim, lai = lai_median, npp = npp_median)
write.csv(all_covariates, file = '/mnt/research/nasabio/data/fia/imputation_covariates.csv', row.names = FALSE)