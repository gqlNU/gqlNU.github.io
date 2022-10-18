rm(list=ls())

library(RColorBrewer)
workdir <- '/Volumes/WorkSpace/OnGitHub/sptmbook/inst/extdata/'
setwd(workdir)
load('copd_simdata.RData')

r <- simdata$y/simdata$e
nareas <- nrow(r)
rho <- rho.sd <- exceed <- rep(0,nareas)
for (i in 1:nareas) {
    m <- arima(r[i,],order=c(1,0,0))
    rho[i] <- m$coef[1]
    rho.sd[i] <- sqrt(m$var.coef[1,1])
    if (rho[i] > 2*rho.sd[i]) exceed[i] <- 1
}

########################################################
####    load and simplify the shapefile
########################################################
require(maptools)
shp <- readShapePoly('England_LAD_2001_simplified.shp')

########################################################
####    plot
########################################################
cutpoints <- c(-1,-0.6,-0.2,0.2,0.6,1)
cutpoints <- quantile(rho)
lv <- cut(rho,cutpoints,labels=1:(length(cutpoints)-1),include.lowest=TRUE)
shadings <- grey(c(0.9,0.75,0.45,0.2))
jpeg(file='ar1.jpg',width=800,height=800,quality=100)
print(plot(shp,col=shadings[lv],lwd=0.8))
intervals <- c('[-0.87, -0.25]', '(-0.25, 0.04]', '(0.04, 0.35]', '(0.35, 0.84]')
legend('topleft',legend=intervals,fill=shadings,cex=2.6,bty='n')
dev.off()
