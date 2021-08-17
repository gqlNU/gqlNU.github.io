###  clear R environment
rm(list=ls())

###  setting working directory 
###  (where data and shapefiles are stored)
workdir <- 'C:/sptmbook/ch8'
setwd(workdir)

###  load libraries
library(R2WinBUGS) #  for running WinBUGS through R
library(maptools)  #  for dealing with shape file
library(spdep)     #  for constructing W

###  load income data
load('income_for_ch8.RData')

###  read Newcastle MSOA shapefile
newcastle <- readShapePoly('Newcastle_MSOAs.shp')
plot(newcastle) #  just checking if everything works

###  construct the W matrix (rook's move)
###  (see chapter 4 for more detail)
nb <- poly2nb(newcastle,queen=FALSE)
spadj <- nb2WB(nb)

###  model file
model.file <- paste0(workdir,'/fig8_1.txt')

###  put all data together for running WinBUGS
wb.data <- list(nhhs=household.income$n,
                nmsoas=household.income$nmsoas,
                y=household.income$income,
                msoa=household.income$msoa)
wb.data <- c(wb.data,spadj)

###  initial values for chain 1
inits1 <- list(alpha=200,sigma.y=50,sigma.S=20)
inits1$S <- rep(0,wb.data$nmsoas)
inits1$S[1] <- 1
inits1$S[2] <- -1
###  initial values for chain 2
inits2 <- list(alpha=700,sigma.y=80,sigma.S=50)
inits2$S <- rep(0,wb.data$nmsoas)
inits2$S[1] <- -1
inits2$S[2] <- 1
inits <- list(inits1,inits2)

###  set parameters to monitor
par2save <- c('alpha','theta','sigma.S','sigma.y','vpc')

###  running WinBUGS
fits <- bugs(data=wb.data,inits=inits,
             parameters.to.save=par2save,
             model.file=model.file,
             n.chains=2,
             n.iter=10000, n.burnin=3000, n.thin=1,
             debug=TRUE, # keep WinBUGS open after running
             # below specifies where WinBUGS is 
             # installed on your machine
             bugs.directory='C:/WinBUGS14',
)

###  display summary of some parameters on screen
fits$summary[c('alpha','sigma.y','sigma.S','vpc'),]

###  save results
save.file <- 'fig8_1_fits.RData'
save(file=save.file,fits,wb.data)