x <- dir()
allfileids <- paste0(rep(c(1000,2000,3000,4000,5000),each=50), '_', 1:50, '.csv')

isthere <- sapply(allfileids, function(i) any(grepl(i, x)))
table(isthere)

paste(as.numeric(which(!isthere)), collapse = ',')
