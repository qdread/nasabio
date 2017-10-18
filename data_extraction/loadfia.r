# Load FIA coordinates
library(dplyr)
fiapnw <- read.csv('/mnt/home/qdr/data/pnw.csv', stringsAsFactors = FALSE)
fiacoords <- fiapnw %>%
	filter(!is.na(ACTUAL_LAT)) %>%
	rename(PLT_CN = CN, lat = ACTUAL_LAT, lon = ACTUAL_LON)