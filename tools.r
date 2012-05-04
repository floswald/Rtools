
# R utilities
# ============





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


# linear map from x \subset [down,up] to [new.down,new.up]
# -----------------------------------------

linear.map <- function(x,down,up,new.down,new.up,plotit=FALSE) { 
	# maps x \subset [down,up] into [new.down,new.up]. new.down > 0!
	stopifnot(x >= down & x <= up)
	stopifnot(new.down>=0 & new.down<new.up)
	rval <- (x - down + new.down)/(up - down)*new.up
	if (plotit) plot(x=x,y=rval,main=paste("linear mapping from [",paste(range(x),collapse=","),"] into [",paste(c(new.down,new.up),collapse=","),"]",sep=""))
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
