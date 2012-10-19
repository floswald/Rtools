# discrete choice with iid uncertainty
# 2 period example
# aim: plot period 1 value function and objective function
# when there is a "kink"
rm(list=ls(all=T))
# questions:
# how does quality of approximation react to changes in number of points used and degree of spline chosen?

library(data.table)
library(splines)
library(ggplot2)
library(gridExtra)
library(reshape)


# cases: sample size n.a and degree of spline degree

low.a  <- -2
high.a <- 2
n.a    <- c(10,25,100)
n.exp  <- length(n.a)
a      <- lapply(1:n.exp,function(i) seq(low.a,high.a,le=n.a[i]))
a      <- lapply(1:n.exp,function(i) a[[i]] - a[[i]][floor(n.a[i]/3)])
low.a  <- lapply(1:n.exp,function(i) head(a[[i]],1))
high.a <- lapply(1:n.exp,function(i) tail(a[[i]],1))
low.y  <- 0
high.y <- 3
n.y    <- 9
y      <- seq(low.y,high.y,le=n.y)
pi.y   <- rep(1/n.y,n.y)
degree <- c(1,2,3)


states <- lapply(1:n.exp,function(i) data.table(expand.grid(a=a[[i]],y=y)))
states <- lapply(1:n.exp,function(i) states[[i]][ ,cash:=a+y])

# period 2: u(c) = -(c - c*)^2
bliss <- 4/5*max(states[[1]][,cash])
ufun <- function(x,bliss=bliss) -(x-bliss)^2 

Bkvalue  <- ufun(x=y,bliss=bliss)
NBkvalue <- lapply(1:n.exp,function(i) ufun(x=states[[i]][,cash],bliss=bliss))
vmax     <- lapply(1:n.exp,function(i) apply(cbind(rep(Bkvalue,each=n.a[[i]]),NBkvalue[[i]]),1,max))
vmax     <- lapply(1:n.exp,function(i) array(vmax[[i]],dim=c(n.a[[i]],n.y)))
evmax    <- lapply(1:n.exp, function(i) vmax[[i]] %*% pi.y)
# matplot(a,vmax[[1]])
# matlines(a,evmax[[1]])

# spline approx
# knot placement function deBoor's averageing rule
knot.select <- function(degree,grid){
# returns a knotvector for a grid of data sites and a spline degree
    n <- length(grid)
    # if (n<(degree+1)*2+1) stop("need at least 2*(degree+1) +1 grid points")
    p <- n+degree+1     # number of nodes required for solving Ax=b exactly
    knots <- rep(0,p)
    knots <- replace(knots,1:(degree+1),rep(grid[1],degree+1))
    knots <- replace(knots,(p-degree):p,rep(tail(grid,1),degree+1))
# this put RcppSimpleGetvalplicity of first and last knot in order
# if there are anough gridpoints, compute intermediate ones
    if (n<(degree+1)) stop("to few grid points for clamped curve")
    if (n>(degree+1)){
        for (j in 2:(n-degree)) knots[j+degree] <- mean(grid[j:(j+degree-1)])
    }
    return(knots)
}

# knots
knotlist        <- lapply(1:n.exp,function(i) lapply(1:3, function(j) knot.select(degree=degree[j],grid=a[[i]]) ) )
names(knotlist) <- paste("points",n.a,sep="")
for (i in 1:n.exp) names(knotlist[[i]]) <- paste("degree",degree,sep="")

# show knot placements 
# produce plots that show knot spaceing
# find multiple knots
count <- lapply(1:n.exp,function(i) lapply(1:3, function(j) as.numeric(table(knotlist[[i]][[j]] ))  ))
# prepare for plotting
count.out <- vector("list",n.exp)
for (i in 1:n.exp) count.out[[i]] <- vector("list",3)
for (i in (1:n.exp)){
	for (j in 1:3){
		for (k in 1:length(count[[i]][[j]])) count.out[[i]][[j]] <- c(count.out[[i]][[j]],1:count[[i]][[j]][[k]])
	}
}

dfs <- lapply(1:n.exp,function(i) lapply(1:3, function(j) data.frame(knots=knotlist[[i]][[j]],mult=count.out[[i]][[j]] ) ) )
ddf <- lapply(1:n.exp,function(i) rbind(dfs[[i]][[1]],dfs[[i]][[2]],dfs[[i]][[3]]))
for (i in 1:length(ddf)) ddf[[i]]$group <- factor(rep(c("deg1","deg2","deg3"),c(nrow(dfs[[i]][[1]]),nrow(dfs[[i]][[2]]),nrow(dfs[[i]][[3]]))))
p1 <- ggplot(data=ddf[[1]],aes(x=knots,y=mult)) + geom_point() + facet_grid(group~.) + opts(title=paste("num data points:",n.a[1]))
p2 <- ggplot(data=ddf[[2]],aes(x=knots,y=mult)) + geom_point() + facet_grid(group~.) + opts(title=paste("num data points:",n.a[2]))
p3 <- ggplot(data=ddf[[3]],aes(x=knots,y=mult)) + geom_point() + facet_grid(group~.) + opts(title=paste("num data points:",n.a[3]))
# grid.arrange(p1,p2,p3,ncol=3)
# all in one data.frame
dddf <- rbind(ddf[[1]],ddf[[2]],ddf[[3]])
dddf$points <- factor(rep(paste(n.a,"points"),c(nrow(ddf[[1]]),nrow(ddf[[2]]),nrow(ddf[[3]]))))
pdf(file="knots.pdf")
print(ggplot(data=dddf,aes(x=knots,y=mult)) + geom_point(size=1) + facet_grid(group~points,scale="free") )
dev.off()



# basis matrix
base        <- lapply(1:n.exp, function(i) lapply(1:3, function(j) splineDesign(x=a[[i]],knots=knotlist[[i]][[j]],ord=degree[j]+1) ) )
names(base) <- names(knotlist)
for (i in 1:n.exp) names(base[[i]]) <- paste("degree",degree,sep="")

# coefficients
b <- lapply(1:n.exp, function(i) lapply(1:3, function(j) solve(base[[i]][[j]], evmax[[i]] ) ) )

# predict values
n.pred <- 110
preda <- lapply(1:n.exp, function(i) seq(low.a[[i]],high.a[[i]],le=n.pred))
pred.values <- lapply(1:n.exp, function(i) lapply(1:3, function(j) splineDesign(x=preda[[i]], knots=knotlist[[i]][[j]], ord=degree[j]+1 ) %*% b[[i]][[j]] ) )
pred.df     <- lapply(1:n.exp, function(i) data.frame(assets=preda[[i]],matrix(unlist(pred.values[[i]]),nrow=n.pred)))
# melt data
mdf <- list()
for (i in 1:n.exp) {
	names(pred.df[[i]]) <- c("assets","deg1","deg2","deg3")
	mdf[[i]] <- melt(pred.df[[i]],id.vars="assets")
	mdf[[i]]$npoints <- factor(paste(n.a[i],"points",sep=""))
}

mdf <- rbind(mdf[[1]],mdf[[2]],mdf[[3]])
pdf(file="kinks.pdf",height=10,width=10)
print( ggplot(data=mdf,aes(x=assets,y=value)) + geom_path()+ facet_grid(variable~npoints))
dev.off()








