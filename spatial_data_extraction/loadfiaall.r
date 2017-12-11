# Load FIA coordinates (entire continental USA)
library(dplyr)
username <- Sys.getenv('USER')
fiaall <- read.csv(paste0('/mnt/home/', username, '/data/allfia.csv'), stringsAsFactors = FALSE)
fiacoords <- fiaall %>%
	filter(!is.na(ACTUAL_LAT)) %>%
	rename(PLT_CN = CN, lat = ACTUAL_LAT, lon = ACTUAL_LON)
