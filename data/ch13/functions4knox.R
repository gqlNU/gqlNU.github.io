create.individual.records <- function(mat,coords,nm=NULL) {
    y <- NULL
    ntimes <- ncol(mat)
    nareas <- nrow(mat)
    if (is.null(nm)) nm <- as.character(1:nareas)
    for (tt in 1:ntimes) {
        for (i in 1:nareas) {
            obs <- mat[i,tt]
            if (obs>0) {
                tmp <- data.frame(rep(tt,obs),rep(coords[i,1],obs),rep(coords[i,2],obs),rep(nm[i],obs))
                y <- rbind(y,tmp)
            }
        }
    }
    return(y)
}


neighbours <- function(dt,ds,eps.t,eps.s,ns=NULL,n2s=NULL) {
    #  when used in the permutation test since the spatial locations are kept fixed
    #  no addition calculation is needed for space
    #  so just enter the values of ns and n2s from the original data and the funciton only works our the temporal neighbours
    if (is.null(ns)) {
    # spatial
        m <- as.matrix(ds)
        n <- nrow(m)
        bs <- rep(0,n)
        for (i in 1:n) bs[i] <- length(which(m[,i]<=eps.s))-1
        ns <- sum(bs)/2
        n2s <- sum(bs^2)/2 - ns
    }
    # temporal
    m <- as.matrix(dt)
    n <- nrow(m)
    bt <- rep(0,n)
    for (i in 1:n) bt[i] <- length(which(m[,i]<=eps.t))-1
    nt <- sum(bt)/2
    n2t <- sum(bt^2)/2 - nt
    out <- c(ns,n2s,nt,n2t)
    names(out) <- c('ns','n2s','nt','n2t')
    return(out)
}

var.knox <-function(nb,n) {
    #  on input, nb is the output from the neighbours function
    #            n is the number of events
    ns <- nb['ns']
    n2s <- nb['n2s']
    nt <- nb['nt']
    n2t <- nb['n2t']
    N <- n*(n-1)/2
    v <- (ns*nt)/N + 4*n2s*n2t/(n*(n-1)*(n-2)) + 4*(ns*(ns-1)-n2s)*(nt*(nt-1)-n2t)/(n*(n-1)*(n-2)*(n-3)) - (ns*nt/N)^2
    return(v)
}

dataframe2matrix <- function(d) {
    mat <- matrix(0,ncol=length(d),nrow=length(d[[1]]))
    for (i in 1:length(d)) mat[,i] <- d[[i]]
    return(mat)
}


