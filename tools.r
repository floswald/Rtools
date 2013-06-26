
# R utilities
# ============

# clustered standard errors
# -------------------------

clx <- 
function(fm, dfcw, cluster){
         # R-codes (www.r-project.org) for computing
         # clustered-standard errors. Mahmood Arai, Jan 26, 2008.
	 
	 # The arguments of the function are:
         # fitted model, cluster1 and cluster2
         # You need to install libraries `sandwich' and `lmtest'
	 
         # reweighting the var-cov matrix for the within model
         library(sandwich);library(lmtest)
         M <- length(unique(cluster))   
         N <- length(cluster)           
         K <- fm$rank                        
         dfc <- (M/(M-1))*((N-1)/(N-K))  
         uj  <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
         vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)*dfcw
         coeftest(fm, vcovCL) }

mclx <- 
function(fm, dfcw, cluster1, cluster2){
         # R-codes (www.r-project.org) for computing multi-way 
         # clustered-standard errors. Mahmood Arai, Jan 26, 2008. 
         # See: Thompson (2006), Cameron, Gelbach and Miller (2006)
         # and Petersen (2006).
	 # reweighting the var-cov matrix for the within model
         
         # The arguments of the function are:
         # fitted model, cluster1 and cluster2
         # You need to install libraries `sandwich' and `lmtest'
         
         library(sandwich);library(lmtest)
         cluster12 = paste(cluster1,cluster2, sep="")
         M1  <- length(unique(cluster1))
         M2  <- length(unique(cluster2))   
         M12 <- length(unique(cluster12))
         N   <- length(cluster1)          
         K   <- fm$rank             
         dfc1  <- (M1/(M1-1))*((N-1)/(N-K))  
         dfc2  <- (M2/(M2-1))*((N-1)/(N-K))  
         dfc12 <- (M12/(M12-1))*((N-1)/(N-K))  
         u1j   <- apply(estfun(fm), 2, function(x) tapply(x, cluster1,  sum)) 
         u2j   <- apply(estfun(fm), 2, function(x) tapply(x, cluster2,  sum)) 
         u12j  <- apply(estfun(fm), 2, function(x) tapply(x, cluster12, sum)) 
         vc1   <-  dfc1*sandwich(fm, meat=crossprod(u1j)/N )
         vc2   <-  dfc2*sandwich(fm, meat=crossprod(u2j)/N )
         vc12  <- dfc12*sandwich(fm, meat=crossprod(u12j)/N)
         vcovMCL <- (vc1 + vc2 - vc12)*dfcw
         coeftest(fm, vcovMCL)}




# grid creator 
# -------------------------

scale.grid <- function(rule="tolower",ub,lb,n.points,breaks=NULL,n.in.breaks=NULL,plotit=FALSE){
# places points in an interval according to rule
	# 2 types of grids: 
	# 1) one scaling for entire interval: require (lb,ub,n)
	# 2) different scaling for different intervals: require (lb,ub,breaks=c(),n.in.breaks=c())
	rules <- c("tolower","ttolower","sinh","toupper","tobreak","ttobreak","tobreak-lin","sinhmid","sinhmid2")
	if (!(rule %in% rules)) stop(paste("your rule",rule,"does not match any of",paste(rules,collapse=","),sep=" "))
	n   <- n.points
	stopifnot(ub > lb)
	if (!is.null(breaks)) stopifnot(lb < head(breaks,1) & ub > tail(breaks,1))
	if (is.null(breaks)){
		# type 1 grid scaling: transform entire grid with same rule
		if (rule=="tolower"){
			# concentrate points at lower bound
			out <- rep(0,n)
			off <- 1
			if (lb<0) off <- 1 - lb #  adjust in case of neg bound
			out[1] <- log(lb + off)
			out[n] <- log(ub + off)
			out    <- seq(from=out[1],to=out[n],le=n)
			out    <- exp( out) - off
			if (plotit) {
				oldpar <- par()$mar
				par(mar=c(15,4,15,2))
				plot(x=out,y=rep(1,n),yaxt="n",ylab="",pch=3,xlab="point allocation",main="exp bunch at lower bound")
				par(mar=oldpar)
			}
			return( out )
		} else if (rule=="ttolower"){
			# concentrate MORE at lower bound
			out <- rep(0,n)
			off <- 1
			if (lb<0) off <- 1 - lb #  adjust in case of neg bound
			out[1] <- log( log(lb + off) + off )
			out[n] <- log( log(ub + off) + off )
			out    <- seq(from=out[1],to=out[n],le=n)
			out    <- exp( exp(out) - off ) - off
			if (plotit) {
				oldpar <- par()$mar
				par(mar=c(15,4,15,2))
				plot(x=out,y=rep(1,n),yaxt="n",ylab="",pch=3,xlab="point allocation",main="asset grid. rule: double exp bunching at lower bound")
				par(mar=oldpar)
			}
			return( out )
		} else if (rule=="sinh"){
			# concentrate around middle of interval
			# inverse sine transform
			out <- sinh(seq(asinh(lb),asinh(ub),le=n))
			if (plotit) {
				oldpar <- par()$mar
				par(mar=c(15,4,15,2))
				plot(x=out,y=rep(1,n),yaxt="n",ylab="",xlab="point allocation",pch=3,main="sinh transform entire grid. contracts if lb < 0 ub")
				par(mar=oldpar)
			}
			return( out )
		} else if (rule=="toupper"){
			# concentrate at upper bound
			out <- rep(0,n)
			out[1] <- exp( lb )
			out[n] <- exp( ub )
			out    <- seq(from=out[1],to=out[n],le=n)
			out <- log(out)
			if (plotit) {
				oldpar <- par()$mar
				par(mar=c(15,4,15,2))
				plot(x=out,y=rep(1,n),yaxt="n",ylab="",pch=3,xlab="point allocation",main="asset grid. rule: concentrate at upper bound")
				par(mar=oldpar)
			}
			return( out )
		}
	} else if (length(breaks)==1) {
		if (rule=="tobreak"){
			# concentrate points towards a unique break
			stopifnot(length(n.in.breaks)==length(breaks)+1)
			stopifnot(n==sum(n.in.breaks))
			below <- scale.grid(rule="toupper",ub=breaks[1],lb=lb,n.points=n.in.breaks[1])
			above <- scale.grid(rule="tolower",ub=ub,lb=breaks[1],n.points=n.in.breaks[2]+1)
			out <- c(below,above[-1])
			if (plotit) {
				oldpar <- par()$mar
				par(mar=c(15,4,15,2))
				plot(x=out,y=rep(1,n),yaxt="n",ylab="",pch=3,xlab="point allocation",main="concentrate towards break")
				par(mar=oldpar)
			}
			return( out )
		} else if (rule=="ttobreak"){
			# concentrate points more towards a unique break 
			stopifnot(length(n.in.breaks)==length(breaks)+1)
			stopifnot(n==sum(n.in.breaks))
			below <- scale.grid(rule="toupper",ub=breaks[1],lb=lb,n.points=n.in.breaks[1])
			above <- scale.grid(rule="ttolower",ub=ub,lb=breaks[1],n.points=n.in.breaks[2]+1)
			out <- c(below,above[-1])
			if (plotit) {
				oldpar <- par()$mar
				par(mar=c(15,4,15,2))
				plot(x=out,y=rep(1,n),yaxt="n",ylab="",pch=3,xlab="point allocation",main="concentrate more towards break")
				par(mar=oldpar)
			}
			return( out )
		} else if (rule=="tobreak-lin"){
			# concentrate points towards a unique break: linear below, exp above
			stopifnot(length(n.in.breaks)==length(breaks)+1)
			stopifnot(n==sum(n.in.breaks))
			below <- seq(from=lb,to=breaks[1],le=n.in.breaks[1])
			above <- scale.grid(rule="ttolower",ub=ub,lb=breaks[1],n.points=n.in.breaks[2]+1)
			out <- c(below,above[-1])
			if (plotit) {
				oldpar <- par()$mar
				par(mar=c(15,4,15,2))
				plot(x=out,y=rep(1,n),yaxt="n",ylab="",pch=3,xlab="point allocation",main="concentrate towards break from above, linear below")
				par(mar=oldpar)
			}
			return( out )
		}
	} else if (length(breaks > 1)) {
		if (rule=="sinhmid") {
			# have an interval in the middle with sinh transform, below and above contracting towards that interval
			stopifnot(length(n.in.breaks)==length(breaks)+1)
			stopifnot(n==sum(n.in.breaks))
			# inverse sine transform in [-1,1]
			low <- scale.grid(rule="toupper",lb=lb,ub=breaks[1],n.points=n.in.breaks[1])
			mid <- sinh(seq(asinh(breaks[1]),asinh(breaks[2]),le=n.in.breaks[2]+1))
			high <- scale.grid(rule="tolower",lb=breaks[2],ub=ub,n.points=n.in.breaks[3]+1)
			out <- c(low,mid[-1],high[-1])
			if (plotit) {
				oldpar <- par()$mar
				par(mar=c(15,4,15,2))
				plot(x=out,y=rep(1,n),yaxt="n",ylab="",xlab="point allocation",pch=3,main=paste("asset grid. rule: sinh transform in [",paste(breaks,collapse=","),"]",sep=""))
				par(mar=oldpar)
			}
			return( out )
		} else if (rule=="sinhmid2") {
			# have an interval in the middle with sinh transform, below and above contracting towards that interval, stronger above
			stopifnot(length(n.in.breaks)==length(breaks)+1)
			stopifnot(n==sum(n.in.breaks))
			# inverse sine transform in [-1,1]
			low <- scale.grid(rule="toupper",lb=lb,ub=breaks[1],n.points=n.in.breaks[1])
			mid <- sinh(seq(asinh(breaks[1]),asinh(breaks[2]),le=n.in.breaks[2]+1))
			high <- scale.grid(rule="ttolower",lb=breaks[2],ub=ub,n.points=n.in.breaks[3]+1)
			out <- c(low,mid[-1],high[-1])
			if (plotit) {
				oldpar <- par()$mar
				par(mar=c(15,4,15,2))
				plot(x=out,y=rep(1,n),yaxt="n",ylab="",xlab="point allocation",pch=3,main=paste("asset grid. rule: sinh transform in [",paste(breaks,collapse=","),"], strong above",sep=""))
				par(mar=oldpar)
			}
			return( out )
		} else {
			stop("scale.grid::you landed in a place with no option")
		}
	}
}


# linear map from x \subset [down,up] to [0,new.up]
# -----------------------------------------

linear.map <- function(x,down,up,new.up,plotit=FALSE) { 
	# maps x \subset [down,up] into [0,new.up]. 
	stopifnot(x >= down & x <= up)
	rval <- (x - down)/(up - down)*new.up
	if (plotit) plot(x=x,y=rval,main=paste("linear mapping from [",paste(range(x),collapse=","),"] into [0",new.up,"]",sep=""))
	return(rval)
}


# nonlinear maps from range(z) to [low,high]
# ------------------------------------------

nonlinear.map <- function(z,low,high,type,plotit=FALSE){
    if (type==1) rval <- ((low + (z-z[1] )/(1 + (z-z[1])/high)))
    if (type==2) rval <- ((low + (z-z[1])^0.6 )/(high + (z-z[1])^0.6))
    if (type==3) rval <- (high*(low + exp(z))/(high + exp(z)))
    if (type==4) rval <- ((low + exp(z))/(1 + (exp(z)/high)))
	if (plotit) {
		plot(x=z,y=rval,main=paste("nonlinear map type: ",type))
		abline(h=np$high.y)
		abline(h=np$low.y)
	}
	return(rval)
}



# get normal copula
# ------------------------------------------

getNormCop <- function(rho,n,Qn= seq(1/n,1-1/n,l=n),cond=FALSE) {

 require(copula)
 cop = normalCopula(param=rho)

 #create the grid
 vals = expand.grid(p=Qn,p2=Qn)

 # apply the copula
 vals$v = dcopula(cop,as.matrix(vals))
 G = array(vals$v,dim=c(n,n))

 # making it conditional
 if (cond) {
   G = t(apply(G,1,function(v) { return(v/sum(v)) }))
 }
 return(G)
}




repmat = function(X,m,n){
##R equivalent of repmat (matlab)
mx = dim(X)[1]
nx = dim(X)[2]
matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}


# spread the array A[1,2,3] to 
# what dim specifies in the order given
# spread(A , c(2) , 10) will insert a dimension
# in the second index of size 10
spread <- function(A, loc, dims) {

  if (!(is.array(A))) {
    A = array(A,dim=c(length(A)))
  }

  adims = dim(A)
  
  # check dimensions
  l = length(loc)
  if (max(loc)> length(dim(A))+l) { 
    stop('incorrect dimensions in spread')
  }

  # total dimension not in order
  sdim = c(dim(A),dims)

  # permutation of dims
  edim = c()
  oi =1        # original dims
  ni =length(dim(A))+1 # new dims
  for (i in c(1: (length(dim(A))+l))) {
    if (i %in% loc) {
      edim = c(edim,ni)
      ni = ni+1
    } else {
      edim = c(edim,oi)
      oi = oi +1
    }
  }

  return( aperm( array(A,dim=sdim),edim)) 
}


# kron.prod 
# =========

# avoiding the curse of dimensionality when 
# multiplying large kronecker products with a vector

# R-implementation of multi_module.f90 by Lars Nesheim
# original code in shared-lib/people/lars/fortran/tools
# based on de Boor (1979) 

# in matrices: fastest index first

kron.prod <- function(y,matrices){
	# INPUT
	# matrices: list(matrix1, matrix2,..., matrixd)
	# 			where dim(matrixj) = c(nj,nj) [SQUARE!]
	# y: value vector 
	# 	 with length(y) = n1*...*nd
	#
	# OUTPUT
	# y1: kron(matrixd,kron(...,(kron(matrix2,matrix1))...))*y
	#
	# check user input
	stopifnot(is.list(matrices))
	# stop if not all of them are square
	allsquare <- lapply(matrices,function(x){dim(x)[1]==dim(x)[2]})
	stopifnot(all(unlist(allsquare)))
	# get sizes
	nmatrices <- length(matrices)
	nall <- length(y)
	nullvec   <- rep(0,nall)
	# compute product for first matrix
	y0    <- y
	stemp <- matrices[[1]]
	n     <- dim(stemp)[1]
	m     <- nall/n
	y1    <- nullvec
	for (i in 1:m){
		y1[m*(0:(n-1)) + i]  <- stemp %*% y0[(n*(i-1)) + (1:n)]
	}
	if (nmatrices > 1){
		# for all other matrices
		for(imat in 2:nmatrices){
			y0    <- y1
			stemp <- matrices[[imat]]
			n     <- dim(stemp)[1]
			m     <- nall/n
			y1    <- nullvec
			for (i in 1:m){
				y1[m*(0:(n-1)) + i]  <- stemp %*% y0[(n*(i-1)) + (1:n)]
			}
		}
	}
	return(y1)
}


# test that in tests/
# =========



# rouwenhorst discretization for AR1
# ----------------------------------
# http://www.karenkopecky.net/rouwenhorst.m
# mu and sigma: mean and standard deviation of error term
# rho: 1st order autocorrelation
# n: number of approximation points
rouwenhorst <- function(rho,sigma,mu,n){
	stopifnot(n>1)
	qu <- (rho+1)/2
	nu <- ((n-1)/(1-rho^2))^(1/2) * sigma
	P  <- matrix(c(qu,1-qu,1-qu,qu),nrow=2,ncol=2)
	if (n>2){
		for (i in 2:(n-1)){
			zeros <- rep(0,i)
			zzeros <- rep(0,i+1)
			P <- qu * rbind(cbind(P,zeros,deparse.level=0),zzeros,deparse.level=0) + (1-qu) * rbind(cbind(zeros,P,deparse.level=0),zzeros,deparse.level=0) + (1-qu) * rbind(zzeros,cbind(P,zeros,deparse.level=0),deparse.level=0) + qu * rbind(zzeros,cbind(zeros,P,deparse.level=0),deparse.level=0)
			P[2:i, ] <- P[2:i, ]/2
		}
	}
	zgrid <- seq(from=mu/(1-rho)-nu,to=mu/(1-rho)+nu,length=n)
	return(list(Pmat=P,zgrid=zgrid))
}




knot.select.old <- function(degree,grid){
# returns a knotvector for a grid of data sites and a spline degree
	grid <- grid[order(grid)]
    n <- length(grid)
    # if (n<(degree+1)*2+1) stop("need at least 2*(degree+1) +1 grid points")
    p <- n+degree+1     # number of nodes required for solving Ax=b exactly
    knots <- rep(0,p)
    knots <- replace(knots,1:(degree+1),rep(grid[1],degree+1))
    knots <- replace(knots,(p-degree):p,rep(tail(grid,1),degree+1))
# this puts multiplicity of first and last knot in order
# if there are anough gridpoints, compute intermediate ones
    if (n<(degree+1)) stop("to few grid points for clamped curve")
    if (n>(degree+1)){
        for (j in 2:(n-degree)) knots[j+degree] <- mean(grid[j:(j+degree-1)])
    }
    return(knots)
}
	


# obtain stationary distribution of transition matrix by iteration
stationary.dist <- function(G,g,n){
	for (i in 1:n){
		g <- G %*% g
	}
	return(g)
}

# simulate n paths of length ntime each from a markov transition matrix trans and initial distribution init.
sim.markov.paths <- function(n,ntime,trans,init=NULL){
	if (is.null(init))	init <- stationary.dist(G=trans,g=rep(1/dim(trans)[1],times=dim(trans)[1]),n=1000)	
	out       <- matrix(NA,nrow=n,ncol=ntime)
	out[,1]   <- sample(1:length(init),size=n,replace=T,prob=init)
	for (i in (1:n)){
		for (it in (2:ntime)){
			out[i,it] <- sample(1:length(init),size=1,replace=T,prob=trans[out[i,it-1], ])
		}
	}
	return(out)
}


# Regression models utilities
# ===========================


# "keep" or "remove" certain variables from a list of summary.lm (or other model) objects
summary.models <- function(model.list,remove=NULL,keep=NULL){
  co <- lapply(model.list,function(x) coef(summary(x)))
  nm <- lapply(model.list,function(x) summary(x)$call)
  names(co) <- nm
  if (is.null(keep)){
	  co <- lapply(co,function(x) x[substr(rownames(x),start=1,stop=8) != remove,])  # keep only variables not named "remove"
  } else {
	  co <- lapply(co,function(x) x[substr(rownames(x),start=1,stop=8) == keep,])  # keep variables named keep
  }
  co <- lapply(co,add.stars)
  return(co)
}


# get grep(which.var) from a summary.lm object
getCoefs <- function(model,which.var,which.cols=1){
	# extracts coefficients "which" from model
	# by default returns only the estimate (col 1)
	x <- coef(summary(model))
	y <- data.frame(var1=grep(which.var,rownames(x),value=TRUE))
	y <- cbind(y,x[grep(which.var,rownames(x)),which.cols])
	names(y) <- c(which.var,colnames(x)[which.cols])
	rownames(y) <- NULL
	return(y)
}


#### add significance stars to a coef(summary.lm(model)) object 

add.stars <- function(model){
	if (!is.matrix(model)) stop( "model needs to be from coef(summary(model)), i.e. a named matrix")
	if ( !"Pr(>|t|)" %in% colnames(model)) stop( "you need to supply the column named 'Pr(>|t|)'")
	df <- as.data.frame(model)
	old <- names(df)
	names(df) <- c("est","std","tval","pval")
	df$sig <- cut(abs(df$pval),breaks=c(0,0.01,0.05,0.1,Inf),right=FALSE,labels=c("***","**","*",""))
	names(df) <- c(old,"")
	return(df)
}


# Time Series Utilities
# ====================

# aggregate zoo object over years
as.year <- function(x) as.numeric(floor(as.yearmon(x)))

# example
# x <- as.yearmon(2000 + seq(0, 23)/12)
# xx <- zoo(seq_along(x), x)
# aggregate(xx, as.year, mean)



