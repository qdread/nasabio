x <- dir()
allfileids <- paste0(rep(c(1000,2000,3000,4000,5000),each=50), '_', 1:50, '.csv')
allfileids <- paste0(rep(as.character(as.integer(c(50000, 75000, 100000, 150000, 200000))),each=50), '_', 1:50, '.csv')
allfileids <- paste0(rep(as.character(as.integer(c(5000, 7500, 10000, 20000, 50000))),each=50), '_', 1:50, '.csv')

allfileids <- paste0('_', 1:150, '.csv')

x <- dir()
allfileids <- paste0('stats_', 1:250, '.r')

x <- dir()
allfileids <- paste0('tdbeta_', rep(1997:2015, times=13), '_', rep(1:13, each=length(1997:2015)), '.r')

x <- dir()
allfileids <- paste0('route_betabaselga_', rep(as.character(as.integer(c(50000, 75000, 100000, 150000, 200000, 300000))),each=20), '_', 1:20, '.csv')

x <- dir()
allfileids <- paste0('betabaselga_', rep(as.character(as.integer(c(1000, 5000, 7500, 10000, 20000, 50000, 75000, 100000, 150000, 200000, 300000))),each=20), '_', 1:20, '.csv')


isthere <- sapply(allfileids, function(i) any(grepl(i, x)))
table(isthere)

paste(as.numeric(which(!isthere)), collapse = ',')


FIA BETA SLICES (194 complete, 56 to do)
171,172,173,174,178,179,182,183,202,203,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250

BBS big radius BETA SLICES (182 complete, 68 to do)
152,153,154,156,157,158,159,160,161,162,165,166,168,170,187,194,195,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250


BBS GAMMA
17,23,24,25,86,87,111,136,151,152,167,168,183,184,188

BBS ELEVATION new
1,2,3,4,5,6,7,13,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,87,88,89,90,91,92,94,95,96,97,98,99,101,103,105,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,130,131,132,133,134,135,137,138,139,140,141,142,143,144,145,146,147,148,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,172,173,174,175,176,177,178,180,181,182,183,184,185,187,188,189,190,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,225,226,228,229,230,231,232,235,238,239,240,241,242,243,245,246,247,248,249,250

BBS betabaselga
19,28,29,30,33,34,39,40,42,43,44,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60

FIA betabaselga
46,47,48,49,50,59,60,63,69,76,77,78,79,80,81,82,83,84,85,86,87,88,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,121,122,133,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220