library(picante)

# Fix match.comm.dist to suppress output
match.comm.dist.fixed <- function (comm, dis) 
{
  res <- list()
  commtaxa <- colnames(comm)
  if (is.null(commtaxa)) {
    stop("Community data set lacks taxa (column) names, these are required to match distance matrix and community data")
  }
  disclass <- class(dis)
  dis <- as.matrix(dis)
  distaxa <- rownames(dis)
  if (is.null(distaxa)) {
    warning("Distance matrix lacks taxa names, these are required to match community and distance matrix. Data are returned unsorted. Assuming that distance matrix and community data taxa columns are in the same order!")
    if (disclass == "dist") {
      return(list(comm = comm, dist = as.dist(dis)))
    }
    else {
      return(list(comm = comm, dist = dis))
    }
  }
  if (!all(distaxa %in% commtaxa)) {
    #print("Dropping taxa from the distance matrix because they are not present in the community data:")
    #print(setdiff(distaxa, commtaxa))
    dis <- dis[intersect(distaxa, commtaxa), intersect(distaxa, 
                                                       commtaxa)]
    distaxa <- rownames(dis)
  }
  if (any(!(commtaxa %in% distaxa))) {
    #    print("Dropping taxa from the community because they are not present in the distance matrix:")
    #   print(setdiff(commtaxa, distaxa))
    res$comm <- comm[, intersect(commtaxa, distaxa)]
  }
  else {
    res$comm <- comm
  }
  if (disclass == "dist") {
    res$dist <- as.dist(dis[colnames(comm), colnames(comm)])
  }
  else {
    res$dist <- dis[colnames(comm), colnames(comm)]
  }
  return(res)
}
assignInNamespace('match.comm.dist', match.comm.dist.fixed, 'picante')