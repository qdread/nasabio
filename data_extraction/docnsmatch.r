# Check whether the fuzzed and unfuzzed are close enough for fia pnw.

# Load fuzzed

library(dplyr)
fiapnw <- read.csv('/mnt/research/nasabio/data/fia/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
  group_by(STATECD, COUNTYCD, PLT_CN, PLOT) %>%
  summarize(lat = LAT_FUZZSWAP[1],
            lon = LON_FUZZSWAP[1])

# Load unfuzzed

pnw_unf <- read.csv('~/data/pnw.csv', stringsAsFactors = FALSE)

pnw_unf <- pnw_unf %>%
	rename(PLT_CN = CN) %>%
	left_join(fiacoords) 
	
pnw_dists <- pnw_unf %>%
	filter(!is.na(lat), !is.na(lon)) %>%
	rowwise %>%
	summarize(dist = sp::spDistsN1(cbind(lon_nad83, lat_nad83), cbind(lon, lat), longlat = TRUE))

table(pnw_dists > 1) # Almost all of them are less than 1 km except for 800 plots.
# However there is something weird where there are a lot more plots in the new versus the old.