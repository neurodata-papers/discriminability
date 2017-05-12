###one sample test for reliability
source('reliability.R')
reliability_test<-function(dist,ids,nrep=100){
  # test whether discriminability differs from 0.5
  N <- dim(dist)[1]
  if (is.null((N))) {
    stop('Invalid datatype for N')
  }
  tr <- mean(rdf(dist, ids))
  
  nr <- rep(0,nrep)  
  for (i in 1:nrep){
    nr[i] <- mean(rdf(dist, ids[sample(N)]))
  }
  result <- list()
  result$rel <- tr
  result$null <- sort(nr)
  result$pv <- sum(nr>tr)/nrep
  return(result)
}

reliability_pipes_test<-function(dist1,dist2,ids,nrep=100){
  # test two discriminability are the same
  N1 <- dim(dist1)[1]
  N2 <- dim(dist2)[1]
  if (is.null(N1) || is.null(N2) ) {
    stop('Invalid datatype for dist1 or dist2')
  }
  if (N1 != N2) {
    stop('The dimension of dist1 and dist2 do not match')
  }
  disct1 <- matrix(0,N1,2)
  disct2 <- matrix(0,N1,2)
  
  for (i in 1:N1){
    disct1[i,] <- dis_vec(dist1[i,],i,ids)
    disct2[i,] <- dis_vec(dist2[i,],i,ids)
  }
  
  tcount <- sum(disct1[,2])
  tdif <- (sum(disct1[,1] * disct1[,2]) - sum(disct2[,1] * disct2[,2])) / tcount
  ##print(sum(disct1[,1] * disct1[,2])/tcount)
  ##print(sum(disct2[,1] * disct2[,2])/tcount)
  
  ndif <- rep(0,nrep)
  ndisct1 <- matrix(0,N1,2)
  ndisct2 <- matrix(0,N1,2)
  for (i in 1:nrep){
    ind <- rbinom(N1,1,0.5) == 1
    for (j in 1:N1){
      ndisct1[ind,] <- disct1[ind,]
      ndisct1[!ind,] <- disct2[!ind,]
      ndisct2[ind,] <- disct2[ind,]
      ndisct2[!ind,] <- disct1[!ind,]
    }
    ndif[i] <- (sum(ndisct1[,1] * ndisct1[,2]) - sum(ndisct2[,1] * ndisct2[,2])) / tcount
  }
  pvalue <- (sum(ndif > abs(tdif)) + 0.5 * sum(ndif == abs(tdif))) / nrep
  return (pvalue)
  
}

dis_vec<-function(distvec,i,ids){
  N <- length(distvec)
    ind <- which(grepl(ids[i],ids))
    rdf <- c()
    count <- 0
    for (j in ind) {
      if (j != i) {
        di <- distvec
        d <- di[j]
        di[ind] <- Inf
        count <- count + 1
        rdf[count] <- 1 - (sum(di[!is.nan(di)] < d) + 0.5*sum(di[!is.nan(di)] == d)) / (N-length(ind))
      }
    }
    return(c(mean(rdf),count))
}


# reliability_bootstrap<-function(dist1,ids1,ids2,nrep=100,obs1=NULL){
#   # test two discriminability are the same
#   N1 <- dim(dist1)[1]
#   if (is.null(N1) ) {
#     stop('Invalid datatype for dist1 or dist2')
#   }
#   
#   rels <- rep(0,nrep)
#   uids1 <- unique(ids1)
#   uids2 <- unique(ids2)
#   
#   if (length(uids1) < length(uids2) ) {
#     print('ID1 is shorter than ID2. Need to bootstrap subjects. obs1 must be provided.')
#   }
#   
#   
#   tids<-ids1
#   tuids1<-uids1
#   for (i in 1:nrep){
#     ids1<-tids
#     uids1<-tuids1
#     tobs<-obs1
#     if (length(uids1) < length(uids2) ) {
#       for (m in 1:(length(uids2) - length(uids1))){
#         t <- runif(1)
#         twoids <- sample(uids1,2)
#         nobs <- obs1[which(grepl(twoids[1],ids1)),]*t + obs1[which(grepl(twoids[2],ids1)),] * (1-t)
#         tobs <- rbind(tobs,nobs)
#         ids1 <- c(ids1, rep(paste('pseudo',m+10000,sep=""),sum(grepl(twoids[1],ids1))))
#       }
#       uids1 <- unique(ids1)
#       dist1 <- as.matrix(dist(tobs,diag = T, upper = T))
#     }
#     
#     
#     
#     ind <- c()
#     tid <- c()
#     uids1 <- sample(uids1,length(uids1))
#     for (j in 1:length(uids2)){
#       nid <- length(which(grepl(uids2[j],ids2)))
#       jid <- which(grepl(uids1[j],ids1))
#       tid <- c(tid, rep(uids2[j], nid))
#       if (length(jid == nid))      {
#         ind <- c(ind, sample(jid,nid,rep=T))
#       } else {
#         ind <- c(ind, sample(jid,nid))
#       }
#       
#     }
#     
#     if (is.numeric(tid)){tid <-tid+1000000}
#     rels[i] <- mean(rdf(dist1[ind,ind], tid))
#   }
#   
#   return(rels)
#   
# }

reliability_bootstrap<-function(ids1,obs1,nrep=100){
  # test two discriminability are the same
  N <- length(ids1)
  if (is.null(N) ) {
    stop('Invalid datatype for dist1 or dist2')
  }
  
  rels <- rep(0,nrep)
  uids1 <- unique(ids1)
  bobs <- obs1
  pduct <- obs1 %*% t(obs1)
  
  
  
  for (i in 1:nrep){
      bids1 <- c()
      cind <- 1    
      
      ind1v<-c()
      ind2v<-c()
      tv<-c()
      
      for (m in 1:(length(uids1))){
        t <- runif(1)
        nmeas <- length(which(grepl(uids1[m],ids1)))
        nmeas1<-0
        nmeas2<-0
        while(min(nmeas1,nmeas2)<nmeas){
        twoids <- sample(uids1,2)
        nmeas1 <- length(which(grepl(twoids[1],ids1)))
        nmeas2 <- length(which(grepl(twoids[2],ids1)))
        }
        ind1<-sample(which(grepl(twoids[1],ids1)),nmeas)
        ind2<-sample(which(grepl(twoids[2],ids1)),nmeas)
        ind1v<-c(ind1v,ind1)
        ind2v<-c(ind2v,ind2)
        tv <- c(tv,rep(t,nmeas))
        
        nobs <- obs1[ind1,]*t + obs1[ind2,] * (1-t)
        bobs[cind:(cind+nmeas-1),] <- nobs
        cind <- cind + nmeas
        bids1 <- c(bids1, rep(paste('pseudo',m+10000,sep=""),nmeas))
      }

      
    distmat1 <- matrix(0,N,N)
    
    

    if (ncol(obs1) > 1){
    for (j in 1:N){
      for(k in j:N){
        dv <- rep(0,N)
        dv[ind1v[j]] <- dv[ind1v[j]] + tv[j]
        dv[ind2v[j]] <- dv[ind2v[j]] + (1-tv[j])
        dv[ind1v[k]] <- dv[ind1v[k]] - tv[k]
        dv[ind2v[k]] <- dv[ind2v[k]] - (1-tv[k])
        distmat1[j,k] <- t(dv) %*% pduct %*% dv
      }
    }
    distmat1 <- sqrt(distmat1)
    distmat1<-distmat1 + t(distmat1) 
    } else {
      bO <- obs1[ind1v,] * tv +  obs1[ind2v,] *(1-tv)
      distmat1 <- as.matrix(dist(bO,diag = T, upper = T))
      distmat1 <- distmat1+ t(distmat1) 
    }
      
    
    rels[i] <- mean(rdf(distmat1,bids1))
  }
  return(rels)
}


