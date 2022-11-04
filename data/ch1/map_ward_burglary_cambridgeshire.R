rm(list=ls());gc()
workdir <- '/Volumes/WorkSpace/ch1'
setwd(workdir)
library(maptools)

#  read the Cambridgeshire ward shapefile in
cam <- readShapePoly('camburgy.shp')

#  extract the relevant data from the shapefile
cam.data <- cam@data
cam.data$burgRate <- cam.data$BURGLARYRA
O <- matrix(cam.data$DBURGLARY,ncol=1)  #  ward-level burglary counts
pop <- matrix(cam.data$DWELLING,ncol=1) #  number of houses

burgRate <- O/pop*1000  #  calculate the burglary rate per 1000 houses

#  extract the centroids of the polygons (for adding the X and + in)
coord <- coordinates(cam)

##########################################################################
####   choropleth map
##########################################################################
cutpoints <- quantile(burgRate)  # obtain min, Q1, median, Q3 and max

#  assign each burglary rate to one fo the four categories defined by quartiles
rate.level <- cut(burgRate,cutpoints,labels=1:4,include.lowest=TRUE)

#  save figure as pdf
pdf.file <- paste(figdir,'fig_1_1.pdf',sep='')
pdf(file=pdf.file,width=10,height=7)

#  define colour shadings
shadings <- grey(c(1,0.7,0.4,0.2))

#  create the choropleth map
plot(cam,col=shadings[rate.level],main='')

#  add legend
intervals <- c('[14.2, 34.9]', '(34.9, 45.6]', '(45.6, 64.1]', '(64.1, 286.7]')
legend('bottomleft',legend=intervals,fill=shadings,cex=1.4,bty='n')

#  add annotation (see p. 4 in the book for explanation)
text(coord[88,1],coord[88,2],'x',cex=2)
text(coord[99,1]+2000,coord[99,2],'+',cex=2)

dev.off()
