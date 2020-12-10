#nb <- poly2nb(peter,queen=FALSE)
#spW <- nb2mat(nb,style='B')

contruct_typeIV_weights <- function(spW,ntimes) {
####    this function constructs the weights for the space-time interaction of Type IV
####    On input:
####       - spW is an NxN spatial weights matrix with diagonals=0
####       - ntimes denotes the number of time points in the time series
	Dsp <- apply(spW,1,sum)
	Qsp <- diag(Dsp) - spW
	Dtm <- c(1,rep(2,ntimes-2),1)
	Qtm <- diag(Dtm)
	for (tt in 1:ntimes) {
		if (tt==1) Qtm[tt,2] <- -1
		if (tt==ntimes) Qtm[tt,tt-1] <- -1
		if (tt>1 & tt<ntimes) Qtm[tt,c(tt-1,tt+1)] <- -1
	}
	QtypeIV <- kronecker(Qsp,Qtm)
	DtypeIV <- diag(diag(QtypeIV))
	WtypeIV <- DtypeIV - QtypeIV
	adj <- list(num=NULL,weights=NULL,adj=NULL)
	n <- nrow(WtypeIV)
	for (i in 1:n) {
		ids <- which(WtypeIV[i,]!=0)
		adj$num <- c(adj$num,length(ids))
		adj$adj <- c(adj$adj,ids)
		adj$weights <- c(adj$weights,WtypeIV[i,ids])
	}
	return(adj)
}

####   Figure 15.4
# spW <- matrix(c(0,1,1,0,1,0,1,0,1,1,0,1,0,0,1,0),byrow=TRUE,ncol=4)
# contruct_typeIV_weights(spW,5)