x <- dir()
allfileids <- paste0(rep(c(1000,2000,3000,4000,5000),each=50), '_', 1:50, '.csv')
allfileids <- paste0(rep(as.character(as.integer(c(50000, 75000, 100000, 150000, 200000))),each=50), '_', 1:50, '.csv')
allfileids <- paste0(rep(as.character(as.integer(c(5000, 7500, 10000, 20000, 50000))),each=50), '_', 1:50, '.csv')

allfileids <- paste0('_', 1:150, '.csv')

isthere <- sapply(allfileids, function(i) any(grepl(i, x)))
table(isthere)

paste(as.numeric(which(!isthere)), collapse = ',')


FIA BETA SLICES (175 complete, 75 to do)
167,168,170,171,172,173,174,175,176,178,179,180,181,182,183,184,185,186,188,190,192,193,196,197,198,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250

BBS ELEV SLICES (132 complete, 18 to do)
34,39,43,44,45,46,61,62,67,68,72,73,93,94,98,104,114,130