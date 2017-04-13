x <- dir()
allfileids <- paste(rep(c(1000,2000,3000,4000,5000),each=50), 1:50, sep='_')

isthere <- sapply(allfileids, function(i) any(grepl(i, x)))
table(isthere)

paste(as.numeric(which(!isthere)), collapse = ',')
