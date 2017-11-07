neighs <- numeric(25)


for (rad in seq(1000,25000,by=1000)) {
  n_rad <- getNeighbors(fiaalbers, radius = rad)
  neighs[rad/1000] <- mean(sapply(n_rad, nrow))
}

getNeighbors <- function(dat, radius) {
  library(spdep)
  idlist <- dnearneigh(coordinates(dat), 0, radius)
  distlist <- nbdists(idlist, coordinates(dat))
  dflist <- list()
  for (i in 1:length(idlist)) {
    if (any(distlist[[i]] <= radius)) {
      dflist[[i]] <- data.frame(idx = idlist[[i]], dist = distlist[[i]][distlist[[i]] <= radius])
    }
    else {
      dflist[[i]] <- NA
    }
  }
  dflist
}

nhb <- getNeighbors(fiaalbers, radius = 25000)


for (rad in seq(1000,25000,by=1000)) {
  neighs[rad/1000] <- mean(sapply(nhb, function(x) ifelse(class(x)=='data.frame', sum(x$dist <= rad), 0)))
}

library(ggplot2)
dat <- data.frame(radius = seq(1000,25000,by=1000), avg_neighbors = neighs)
ggplot(dat, aes(radius, avg_neighbors)) + geom_point()
