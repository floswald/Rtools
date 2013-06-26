# R howto cheatsheet
# ==================

# Author florian.oswald@gmail.com
# let me know if you find any errors.

# R base plot
# ===========

#### Generate k random walks across time {0, 1, ... , T}. 
#### code from http://www.r-bloggers.com/a-plot-of-250-random-walks/

T <- 100
k <- 250
initial.value <- 10
GetRandomWalk <- function() {
  # Add a standard normal at each step
  initial.value + c(0, cumsum(rnorm(T)))
}
# Matrix of random walks. each column is an independent random walk.
values <- replicate(k, GetRandomWalk())
# Create an empty plot
dev.new(height=7, width=10)
plot(0:T, rep(NA, T + 1), main=sprintf("%s Random Walks", k),
     xlab="time", ylab="value",
     ylim=10 + 4.5 * c(-1, 1) * sqrt(T))
mtext(sprintf("%s%s} with initial value of %s",
              "Across time {0, 1, ... , ", T, initial.value))
for (i in 1:k) {
  lines(0:T, values[ , i], lwd=0.25)
}
for (sign in c(-1, 1)) {
	curve(initial.value + sign * 1.96 * sqrt(x), from=0, to=T,
        n=2*T, col="darkred", lty=2, lwd=1.5, add=TRUE) # use the curve() function
}
legend("topright", "1.96 * sqrt(t)", bty="n", lwd=1.5, lty=2, col="darkred")



# ggplot
# ======

library(ggplot2)	# see http://had.co.nz/ggplot2/

#### make the random walks graph with ggplot.

vdf <- as.data.frame(cbind(1:nrow(values),values))	# make a data.frame from matrix 'values'. attach a first column with a time index to it.
names(vdf) <- c("time",paste("trial",1:k,sep=""))	# fix column names
head(vdf)
library(reshape)	# ggplot works best with molten data. see http://had.co.nz/reshape/
mdf <- melt(vdf,id.vars="time")
head(mdf)
p1 <- ggplot(data=mdf, aes(x=time,y=value,group=variable)) + geom_line(size=0.25)	 # type p1.
# define a function that returns plus/minus the 95% CI of a random walk at time x. this CI is a funciton of time (x).
sefun <- function(x,initial.value,plus=TRUE){ 
	if(plus)  return( initial.value + 1.96 * sqrt(x) ) 
	if(!plus) return( initial.value - 1.96 * sqrt(x) ) }	

# apply the positive and negative confidence bound as a layer to the existing plot.
p1 <- p1 + stat_function(fun = function(x) {sefun(x,initial.value=initial.value)}, color="red", size=0.25)	# stat_function applies fun at each x.
p1 <- p1 + stat_function(fun = function(x) {sefun(x,initial.value=initial.value,plus=FALSE)}, color="red", size=0.25)	# plus=FALSE .



#### plot dates: just have x in as.Date() format. see http://had.co.nz/ggplot2/scale_date.html

df <- data.frame(date = seq(Sys.Date(), len=100, by="1 day")[sample(100, 50)], price = runif(50) )	# randomly choose 50 dates out of 100, and 50 prices. so: ggplot does not need equally spaced dates!
df <- df[order(df$date),]	# bring into ascending order
qplot(data=df, x=date, y=price, geom="line")

#### arrange plots on layout similar to par(mfrow)

library(gridExtra)
dframe <- data.frame(group=rep(factor(LETTERS[1:4]),25))
dframe$x <- rnorm(100)
dframe$y <- rnorm(100)
dframe$z <- rnorm(100)
px <- ggplot(dframe,aes(x,y))
pz <- ggplot(dframe,aes(z,y))
# vertical arrange
p1 <- px + geom_point() + facet_grid(group~.) + coord_equal() + stat_smooth(method="lm")
p2 <- pz + geom_point() + facet_grid(group~.) + coord_equal() + stat_smooth(method="lm")
grid.arrange(p1,p2,ncol=2)
# horizontal arrange
p1 <- px + geom_point() + facet_grid(.~group) + coord_equal() + stat_smooth(method="lm")
p2 <- pz + geom_point() + facet_grid(.~group) + coord_equal() + stat_smooth(method="lm")
grid.arrange(p1,p2,ncol=1)
# can also do grid.arrange(p1,p2,p3,p4,ncol=2,nrow=2) if got 4 different plots to be arranged in panels.


#### fully worked example: arrange several different types of plots into one output. densities. uses viewport functionality.
#### adaptation of http://www.r-bloggers.com/graphics-2012-05-14-220400/
#### have a look at that website. code below has more comments.

library(grid)
data(iris)
x <- jitter(iris[,c('Sepal.Length')])
y <- jitter(iris[,c('Sepal.Width')])
z <- factor(iris[,c('Species')])

# viewport needs a function where we specify row and col position of a plot:
vplayout <- function(x, y) {viewport(layout.pos.row = x, layout.pos.col = y) }
	
df <- data.frame(x,y,z) # put iris data into a dataframe

# this will be the top plot: two density estimates.
p1 <- ggplot(data = df) # base plot. no layers => no plot if you type p1 now.
p1 <- p1 + geom_density(aes(x = x, y = ..count.., fill = z, alpha=0.4)) # add first layer. density (rather "count") split by z. (fill densities with different colors, one color for each z). alpha is an opacity parameter. type p1 now.
p1 + geom_density(aes(x = x, y = -..count.., fill=z),  position = "stack", alpha=0.4)  # add second layer.  plot negative density ("-count") and stack densities on top of each other.  tells you distribution by value of z (upper part) and overall distribution of values of x.
p1 <- p1 + geom_density(aes(x = x, y = -..count.., col=z), fill="#CCCCCCCC", position = "stack", alpha=0.4)  # type p1 now. change the filling color to grey, and the lines to z.
p1 <- p1 + opts(legend.position = "none")	# remove legend. type p1 now.

# second density plot, oriented vertically (hence the 'coord_flip()') at the end
p2 <- ggplot(df) 
p2 <- p2 + geom_density(aes(x = y, y = -..count.., col=z), fill="#CCCCCCCC", position = "stack") + geom_density(aes(x = y, y = ..count.., fill = z, alpha=0.4)) + opts(legend.position = "none") # type p2 now. this is identical to before.
p2 <- p2 + coord_flip()	# now it's flipped.

# finally the main x/y plot 
p3 <- ggplot(df,aes(x = x, y = y, col=z)) + geom_point() + opts(legend.position = c(1.2,1.2)) # type p3!

# now we use viewport to squash those three plots together in a nice picture:
grid.newpage()	# create a new grid.
pushViewport(viewport(layout = grid.layout(5, 5))) # a 5 by 5 grid
print(p1, vp=vplayout(1,1:4))   # we print the p1 to row 1, columns 1:4 of grid
print(p3, vp=vplayout(2:5,1:4)) # print the main graph to region row(2:5), cols(1:4)
print(p2, vp=vplayout(2:5,5))   # print p3 to region row(2:5), col(5)


### arrange ggplot legend
theme(legend.position="top")	# put on top of graph, not on right side.



# read data
# =========

#### from disk
dat <- read.csv("path/to/csv-file.csv",header=TRUE)	# read from a local csv file. see help(read.table)
dat <- read.dta("path/to/dta-file.dta")	# in library(foreign)

#### from html table (only or local)
# scrape data and transform 123.455 into 123455 for several cols of data.frame
library(XML)
myurl       <- "~/Dropbox/git/Rstuff/scrape/hp.html"
tables      <- readHTMLTable(myurl)
mdat        <- tables[[2]]
mdat        <- mdat[-1,]
mdat        <- mdat[,-7]
names(mdat) <- c("Month","Detached","SemiDetached","Terraced","Flat","All")
mdat$Month  <- as.Date(paste("01 ",mdat$Month,sep=""),format="%d %b %Y")
mdat[,-1]   <- apply(mdat[,-1],2,function(x) as.numeric(gsub(",","",as.character(x))))	# get rid of commas as thousand sign
myurl       <- "~/Dropbox/git/Rstuff/scrape/hp.html"
tables      <- readHTMLTable(myurl)



# data.table useage
# =================
library(data.table)
dtab <- data.table(x=rep(1:3,each=3),y=rnorm(9),z=rnorm(9))

# compute summary over many columns
qtab <- dtab[,lapply(.SD,mean,na.rm=TRUE),by="x",.SDcols=2:ncol(dtab)]
# .SD means "sub data.table"


# data.frame usage
# ===============

# drop/keep columns
dframe       <- data.frame(factor1=factor(rep(1:3,each=4)),factor2=factor(rep(1:4,each=3)),values1=rnorm(12),values2=rnorm(12))
head(dframe)
dframe       <- dframe[,!names(dframe) %in% c("factor2","values2")]
sandp        <- sandp[,!(names(sandp) %in% c("CA.Los.Angeles","CA.San.Diego","CA.San.Francisco","FL.Miami","FL.Tampa","Composite.20","YEAR"))]

# cumulative sum over columns of data.frame
cumdf <- apply(dframe,2,cumsum)



# string manipulation
# =========================

# split strings
long   <- c("CA Los.Angeles","FL Tampa", "CT New.Haven", "NY New.York")
splong <- strsplit(long,split=" ")	# split at space
states <- lapply(1:4,function(x) splong[[x]][1])	# states are first element
cities <- lapply(1:4,function(x) splong[[x]][2])	# cities are second
cities <- sub("\\."," ",cities)	# substitute . with white space
# trailing white space
longw   <- c(" CA Los.Angeles "," FL Tampa ", " CT New.Haven ", " NY New.York ")
longend <- sub("\\s+$", '',longw)	# remove white space at end of string
longbeg <- sub("^\\s+", '',longw)	# remove white space at beg of string
# get rid of colon in 'quarter' field
mystr <- "Has:colon"
mystr <- sub('\\:',' ',mystr)
# remove period from beginning of word
gsub("\\.([[:alpha:]])","\\1",".WordWithAPeriod")

# useful tricks
# -------------

# NOT in operator: define 
'%ni%' <- Negate('%in%')



# time series manipulations
# =========================

# TODO: look at library(lubridate)

#### read data
dat   <- read.csv("~/Dropbox/bankruptcy/data/sandp/raw/case-shiller.csv",stringsAsFactors=FALSE)
tsdat <- ts(dat$Composite.10, start = 1987, frequency=12)

#### some plots of a time series: level, first difference, autocorrelation and partial autocorrelation plots
par(mfrow=c(2,2))	# set up a 2 by 2 plotting device
plot(tsdat)
plot(diff(tsdat,differences=1))
acf(tsdat)
pacf(tsdat)
par(mfrow=c(1,1))	# reset the device

#### built-in additive decomposition of time series. uses a moving average by default.
plot(decompose(tsdat))

#### what are first and second differences of the log of this?
par(mfrow=c(2,2))
plot(log(tsdat))
plot(diff(log(tsdat),differences=1))
plot(diff(log(tsdat),differences=2))
par(mfrow=c(1,1))

#### detrend a series manually
plot(tsdat,ylim=c(-50,240),ylab="Case Shiller Index") # plot series in a big enough window
require(zoo)	# load zoo for coredata()
d <- lm(coredata(tsdat) ~ index(tsdat))               # fit a linear time trend
abline(d)                                             # add regression line to the plot
lines(d$residuals ~ index(tsdat),type="l",col="red")  # add the detrended values (i.e. the reg residuals)
abline(h=0,col="red")                                 # add a horizontal reference line at zero
legend(x="topleft",c("observed","detrended"),col=c("black","red"),lty=c(1,1))

#### zoo: make aggregate dates from monthly to quarterly
library(zoo)
mdat <- data.frame(month=seq(Sys.Date(),len=100,by="1 month"),value=rnorm(100)) # notice how the seq() method for dates creates a "Date" object
str(mdat)
head(mdat)
z <- zoo(mdat$value,order.by = mdat$month)                                      # creates a zoo object
qdat <- aggregate(z,as.yearqtr,mean)                                            # aggregate within quarter
head(qdat)

#### zoo: go from quarterly (ad hoc) "20113" to date class "2003-09-01"
library(zoo)
library(data.table) 	# if want to use integer based date class IDate
as.IDate((as.yearqtr("20041",format="%Y%q")))

#### xts: more powerful aggregation possibilities. subsetting by dates is a charm.
require(xts)
# from a daily zoo to weekly xts
tdat <- data.frame(date = seq(Sys.Date(), len=386, by="1 day"), price = rnorm(386)) # daily data
head(tdat)
zdat <- zoo(tdat$price, order.by = tdat$date)
head(zdat)
mdat <- aggregate(zdat, as.yearmon, "mean", na.rm=T)	# use as.yearmon from zoo
head(mdat)
wdat <- apply.weekly(zdat, mean)

# form a xts directly from a data.frame. 
mdat <- data.frame(date=seq(Sys.Date(), len=214, by="1 day"),price=rnorm(214))
mdat
mxts <- as.xts(mdat[,-1],order.by=mdat$date)
mxts
mxts["2013"]                                                     # subsetting by CCYY-MM-DD HH:MM:SS or date characters!
mxts["2012-10-18/2013-01-08"]                                    # dates from / to
first(mxts, "1 week")                                            # first week of data
first(mxts, "3 week")                                            # first 3 weeks of data
first( last(mxts, "2 week"), "3 days")                           # first 3 days of the last 2 weeks
mxts[paste(seq(Sys.Date(),len=386, by="2 day"))]                 # every other day
mxts[paste(seq(Sys.Date(),len=386, by="1 day")[sample(386,15)])] # 15 random days

##### useful xts methods
plot(mxts,major.ticks='months',minor.ticks=FALSE) # plot method
periodicity(mxts)                                 # prints periodicity of data
endpoints(mxts, on = "weeks")                     # prints days where a week ends
endpoints(mxts, on = "months")                    # prints days where a month ends
to.period(mxts, 'months')                         # changes periodicity to monthly, giving first, last and min and max values
apply.weekly(mxts, FUN=mean)                      # apply functions by a certain time frame
apply.monthly(mxts, FUN=max)
apply.quarterly(mxts, FUN=quantile)
apply.yearly(mxts, median)
period.apply(mxts, INDEX=endpoints(mxts, on="quarters"), FUN=max)	# apply functions to user defined periods
period.apply(mxts, INDEX=endpoints(mxts, on="quarters"), FUN=median)




