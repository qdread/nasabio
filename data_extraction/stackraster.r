library(raster)
vars <- c('tmin','tmean','tmax','ppt')
for (v in vars) {
	vstack <- stack(dir(pattern=v))
	writeRaster(vstack, paste('prism',v,'stack_mi.tif', sep='_'), format='GTiff')
}