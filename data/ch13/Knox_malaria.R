# source('/Volumes/WorkSpace/OnGitHub/sptmbook/inst/scripts/knox_malaria.R')

rm(list=ls())
library(rgdal)
library(surveillance)
library(geosphere)  #  to calculate distance between two points in long/lat
source('/Volumes/WorkSpace/OnGitHub/sptmbook/inst/scripts/functions4knox.R')

########################################
####   load shapefile and construct spatial W matrix
########################################
shp <- readOGR('/Volumes/WorkSpace/OnGitHub/sptmbook/inst/extdata/GulMalaria[1].shp')
shpLL <- spTransform(shp,'+init=epsg:4326')

dat <- dataframe2matrix(slot(shpLL,'data')[c(sapply(1:12,function(x){paste('INC',x,sep='')}))])
coords <- coordinates(shpLL)
nm <- as.character(slot(shpLL,'data')$LOCATION)
colnames(coords) <- c('long','lat')
nareas <- nrow(coords)

####   calculate distances across the polygons
fn <- function(position1,position2) distm(c(position1[1],position1[2]),c(position2[1],position2[2]),fun=distHaversine)  #  in metre
d <- proxy::dist(coords,method=fn)/1000
summary(d)
hist(d)

####   space and time thresholds
eps.ss <- c(0,2,6,10,20)
eps.ts <- c(0,1,2,4)
neps.ss <- length(eps.ss)
neps.ts <- length(eps.ts)

y <- create.individual.records(dat,coords,nm)
names(y) <- c('year','long','lat','area')
y$area <- as.character(y$area)

####   calculate distances (in metre)
distance.metre <- proxy::dist(y[,c('long','lat')],method=fn)  #  takes a while to compute

ds <- distance.metre/1000

n <- length(y$year)
####   quantities from the real data
dt <- dist(y$year)
real.test.stats.std <- 
  ns <- 
  n2s <- 
  p.values <- 
  real.test.stats <- 
  mu.real <- 
  v.real <- array(0,c(neps.ss,neps.ts))

for (i in 1:neps.ss) {
    for (j in 1:neps.ts) {
        k <- knox(dt=dt,ds=ds,eps.s=eps.ss[i],eps.t=eps.ts[j],simulate.p.value=FALSE)
        real.test.stats[i,j] <- k$table[1,1]
        mu <- k$null.value
        nb <- neighbours(dt=dt,ds=ds,eps.t=eps.ts[j],eps.s=eps.ss[i])
        ns[i,j] <- nb['ns']
        n2s[i,j] <- nb['n2s']
        v <- var.knox(nb,n=n)
        real.test.stats.std[i,j] <- (k$table[1,1]-mu)/v
    }
}

set.seed(1)
####   permutations
nperms <- 999
test.stats.std <- test.stats <- array(0,c(nperms,neps.ss,neps.ts))
sim.times <- array(0,c(n,nperms))
for (isim in 1:nperms) {
    sim.times[,isim] <- sample(y$year)
    dt.sim <- dist(sim.times[,isim])
    for (i in 1:neps.ss) {
        for (j in 1:neps.ts) {
            k <- knox(dt=dt.sim,ds=ds,eps.s=eps.ss[i],eps.t=eps.ts[j],simulate.p.value=FALSE)
            test.stats[isim,i,j] <- k$table[1,1]
            mu <- k$null.value
            nb.sim <- neighbours(dt=dt.sim,ds=ds,eps.t=eps.ts[j],eps.s=eps.ss[i],ns=ns[i,j],n2s=n2s[i,j])
            v <- var.knox(nb.sim,n=n)
            test.stats.std[,i,j] <- (test.stats[isim,i,j] - mu)/v
        }
    }
    print(isim)
}

for (i in 1:neps.ss) {
    for (j in 1:neps.ts) {
        p.values[i,j] <- (length(which(test.stats[,i,j]>=real.test.stats[i,j]))+1)/(nperms+1)
    }
}
t(p.values)

M <- rep(0,nperms)
for (isim in 1:nperms) {
    M[isim] <- max(test.stats.std[isim,,])
}
(r <- max(real.test.stats.std))
(max(M))
(length(which(M>=r))+1)/(nperms+1)
