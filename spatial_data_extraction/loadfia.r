# Load FIA coordinates
library(dplyr)
username <- Sys.getenv('USER')
fiapnw <- read.csv(paste0('/mnt/home/', username, '/data/pnw.csv'), stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
	filter(!is.na(ACTUAL_LAT)) %>%
	rename(PLT_CN = CN, lat = ACTUAL_LAT, lon = ACTUAL_LON)
