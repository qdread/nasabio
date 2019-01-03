# FIA fuzzed coordinates for Pacific Northwest.

fiapnw <- read.csv('C:/Users/Q/Dropbox/projects/nasabiodiv/finley_trees_pnw_2015evaluations_feb14_2017.csv', stringsAsFactors = FALSE)
fiacoords_fuzzed <- fiapnw %>%
  group_by(PLT_CN, STATECD, COUNTYCD, PLOT) %>%
  summarize(lonfuzz = LON_FUZZSWAP[1],
            latfuzz = LAT_FUZZSWAP[1])
write.csv(fiacoords_fuzzed, 'C:/Users/Q/Dropbox/projects/nasabiodiv/fia_unfuzzed/fia_pnw_coords_fuzzed.csv', row.names = FALSE)
