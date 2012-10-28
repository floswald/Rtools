
# R script accompanying quasi-concave, concave etc note

if(!require(rgl)) install.packages('rgl')
require(rgl)

print("welcome to this R demo. This program will produce a couple of windows, which you can arrange on your screen. Once a window is active (you click inside), you can grag the displayed object with the mouse and rotate it around.")

# setup data
nx <- 50	# points into x direction
ny <- 50	# points into y

x <- seq(0,1,length=50)
y <- seq(0,1,length=50)

df <- expand.grid(x=x,y=y)

# DGPs
cobb.douglas <- function(x,a,b) return( x[1]^a * x[2]^b )

# plotting data
cd1 <- apply(df,1,cobb.douglas,a=0.25,b=0.25)
cd2 <- apply(df,1,cobb.douglas,a=2,b=2)

# put in a matrix
cd.heights1 <- matrix(cd1,nrow=nx,ncol=ny)
cd.heights2 <- matrix(cd2,nrow=nx,ncol=ny)

# concave cobb douglas
open3d()
surface3d(x,y,cd.heights1,color="red",alpha=0.8)
axes3d(labels=FALSE,tick=FALSE)
title3d("concave cobb douglas")

# concave cobb douglas with planes
open3d()
surface3d(x,y,cd.heights1,color="red",alpha=0.8)
axes3d(labels=FALSE,tick=FALSE)
level <- seq(0.25,1.75,by=0.25)
planes3d(a=0,b=0,c=-1,d=level,color="blue",alpha=0.8)
title3d("concave cobb douglas with planes")

contour(x=x,y=y,z=cd.heights1,levels=level,main="contours of concave cobb-douglas")

# quasi-concave cobb douglas
open3d()
surface3d(x,y,cd.heights2,color="red")
axes3d(labels=FALSE,tick=FALSE)
points3d(x=c(1,1),y=c(0,0.8),z=c(cobb.douglas(x=c(1,0),a=2,b=2),cobb.douglas(x=c(1,0.8),a=2,b=2)),size=4)
segments3d(x=c(1,1),y=c(0,0.8),z=c(cobb.douglas(x=c(1,0),a=2,b=2),cobb.douglas(x=c(1,0.8),a=2,b=2)),lwd=2)
title3d("quasi-concave cobb douglas")

open3d()
surface3d(x,y,cd.heights2,color="red")
axes3d(labels=FALSE,tick=FALSE)
level <- seq(0.05,1.95,by=0.25)
planes3d(a=0,b=0,c=-1,d=level,color="blue",alpha=0.8)
title3d("quasi-concave cobb douglas with planes")

dev.new()
contour(x=x,y=y,z=cd.heights2,levels=level,main="contours of quasi-concave cobb-douglas")

print("end of demo")
