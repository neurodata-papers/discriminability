rdf <- function(dist, ids) {
  N <- dim(dist)[1]
  if (is.null((N))) {
    stop('Invalid datatype for N')
  }
  
  uniqids <- unique(as.character(ids))
  countvec <- vector(mode="numeric",length=length(uniqids))
  
  for (i in 1:length(uniqids)) {
    countvec[i] <- sum(grepl(uniqids[i], ids)) # total number of scans for the particular id
  }
  
  scans <- max(countvec) # assume that we will worst case have the most
  rdf <- array(NaN, N*(scans-1)) # initialize empty ra
  
  count <- 1
  for (i in 1:N) {
    ind <- which(grepl(ids[i],ids)) # all the indices that are the same subject, but different scan
    for (j in ind) { 
      if (j != i) { # if j != i, then we want j to have a close distance to i, and estimate where it ranks 
        di <- dist[i,] # get the entire ra for the particular scan
        di[ind] <- Inf # don't want to consider the particular scan itself
        d <- dist[i,j] # the distance between the particular scan of a subject and another scan of the subject
        rdf[count] <- 1 - (sum(di[!is.nan(di)] < d) + 0.5*sum(di[!is.nan(di)] == d)) / (N-length(ind)) # 1 for less than, .5 if equal, then average
        count <-  count + 1
      }
    }
  }
  return(rdf[1:count-1]) # return only the occupied portion
}

mnr <- function(rdf, remove_outliers=TRUE, thresh=0, output=FALSE) {
  if (remove_outliers) {
    mnr <- mean(rdf[which(rdf[!is.nan(rdf)] > thresh)]) # mean of the rdf
    ol <- length(which(rdf<thresh))
    if (output) {
      print(paste('Graphs with reliability <',thresh,'(outliers):', ol))
    }
  } else {
    ol <- 0
    mnr <- mean(rdf[!is.nan(rdf)])
  }
  nopair <- length(rdf[is.nan(rdf)])
  if (output) {
    print(paste('Graphs with unique ids:',nopair))
    print(paste('Graphs available for reliability analysis:', length(rdf)-ol-nopair))
    print(paste('MNR:', mnr))
  }
  return(mnr)
}
