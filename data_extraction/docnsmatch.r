# Check whether the fuzzed and unfuzzed are close enough for fia pnw.

# Load fuzzed

library(dplyr)
fiapnw <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
#fiapnw <- read.csv('/mnt/research/nasabio/data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

# Load unfuzzed

#pnw_unf <- read.csv('~/data/pnw.csv', stringsAsFactors = FALSE)
pnw_unf <- read.csv('~/FIA/pnw.csv', stringsAsFactors = FALSE)

pnw_unf <- pnw_unf %>%
	rename(PLT_CN = CN) %>%
	left_join(fiacoords) 
	
pnw_dists <- pnw_unf %>%
	filter(!is.na(lat), !is.na(lon)) %>%
	rowwise %>%
	summarize(dist = sp::spDistsN1(cbind(lon_nad83, lat_nad83), cbind(lon, lat), longlat = TRUE))

table(pnw_dists > 1) # Almost all of them are less than 1 km except for 800 plots.
# However there is something weird where there are a lot more plots in the new versus the old.

pnw_unf$fuzzdist <- NA
pnw_unf$fuzzdist[!is.na(pnw_unf$lon)] <- as.numeric(pnw_dists$dist)

table(pnw_unf$INVYR[pnw_unf$fuzzdist > 10])

# Figure out the distribution of plots.

with(fiapnw, nrow(unique(cbind(STATECD, COUNTYCD, PLOT))))
with(pnw_unf, nrow(unique(cbind(STATECD, COUNTYCD, PLOT))))

matches <- fiacoords$PLT_CN %in% pnw_unf$PLT_CN

plot(fiacoords$lon, fiacoords$lat, col = as.numeric(matches))
with(fiacoords[matches,], plot(lon, lat))
with(fiacoords[!matches,], plot(lon, lat))

matches2 <- fiapnw$PLT_CN %in% pnw_unf$PLT_CN

# The missing plots are all over the place and are distributed kind of randomly
with(fiapnw[!matches2, ], table(MEASYEAR))

# Make list of CNs of the 2089 missing plots and export as CSV
missingcns <- with(fiacoords[!matches, ], unique(PLT_CN))

write.csv(data.frame(CN = missingcns), file = 'C:/Users/Q/Dropbox/projects/nasabiodiv/missingcns.csv', row.names = FALSE)
