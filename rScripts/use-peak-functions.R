##======================================
##	use-peak-functions.R
##	Sample usage of peak functions
##======================================
setwd("X:/myGit/workingWithYicheng/rScripts/")
# assuming same level of directory, 
# source the function files

source ("../utils/peak-functions.R")



integrate(P, -1, 1, mean=0, sd=1)

# Using the integrate function 
# and probability density function for normal "P"
# The answer is same as "dnorm"
# integrate between 1 standard deviation

integrate(dnorm, -1, 1, mean=0, sd=1)


# integrate(P,-5,5, mean=0, sd=1)   # should be close to 1.0



# Setup test Gaussian with mean 0.0 and standard deviation 1.0
Delta <- 0.01
x <- seq(-5,5, by=Delta)
y <- P(x, 0.0, 1.0)




# peak location:

x[which(peaks(y))]



par(mfrow=c(3,1), oma=c(2,2,2,2))

##### Gaussian
plot(x,y, type="l", lwd=3,
     main="Single Gaussian",
     xlab="x", ylab="y")
abline(h=0, v=1, lty="dotted")

##### 1st Derivative
# According to Abramowitz & Stegun, inflection points are at +/- sigma.
derivative1 <- Deriv1(x,y)

plot(derivative1$x,derivative1$y, type="l",
     main="1st Derivative", xlab="x", ylab="y'")
abline(h=0, v=0, lty="dotted")

peaks.Deriv1   <- peaks( derivative1$y, span=3)
valleys.Deriv1 <- peaks(-derivative1$y, span=3)

points( derivative1$x[peaks.Deriv1], derivative1$y[peaks.Deriv1],
        pch=19, col="red")
points( derivative1$x[valleys.Deriv1], derivative1$y[valleys.Deriv1],
        pch=19, col="blue")

# Approximate location of peak and valley
derivative1$x[peaks.Deriv1]
derivative1$x[valleys.Deriv1]

s.approx <- (derivative1$x[valleys.Deriv1][1] - derivative1$x[peaks.Deriv1][1])/2

#	This is standard deviation!!
s.approx


##============================================
##	Now let's try normal cell, D.I. values
##============================================

y.DI <- P(x, 1.0, 2.0)
y <- y.DI
var.est <- getVarianceViaDeriv (x, y.DI)
mean.est <- getMeanViaDeriv(x, y.DI)

var.est
mean.est

##===============================================
##	Now let's try the mixture of two guassian
##===============================================

Delta <- 0.01
x <- seq(0.0,10.0, by=Delta)
y1 <- P(x, 3.75, 0.75)
y2 <- P(x, 6.00, 0.50)
y <- y1 + y2
plot(x,y)

var.est <- getVarianceViaDeriv (x, y)
var.est
mean.est <- getMeanViaDeriv(x, y)
mean.est


##=================================================
##	Now let's try the mixture of three guassians
##=================================================

Delta <- 0.01
x <- seq(0.0,10.0, by=Delta)
y1 <- P(x, 3.75, 0.75)
y2 <- P(x, 6.00, 0.50)
y3 <- P(x, 8.00, 0.6)
y <- y1 + y2 + y3
plot(x,y)

var.est <- getVarianceViaDeriv (x, y)
var.est
mean.est <- getMeanViaDeriv(x, y)
mean.est


##=================================================
##	Three guassians with different proportion
##=================================================

y <- 0.5*y1 + 0.4*y2 + 0.1*y3
var.est <- getVarianceViaDeriv (x, y)
var.est
mean.est <- getMeanViaDeriv(x, y)
mean.est


##========================================================
##	Three guassians mimicing exfoliative cytology data 
##	Fail when proportions are skewed,
##========================================================

Delta <- 0.01
x <- seq(0.0,10.0, by=Delta)
y1 <- P(x, 1.0, 0.75)
y2 <- P(x, 2.0, 0.50)
y3 <- P(x, 3.6, 0.6)
y <- y1 + y2 + y3

y <- 0.8*y1 + 0.15*y2 + 0.05*y3
plot(x,y)


setwd("X:/myGit/workingWithYicheng/fuyue\ data/CSV/OLK")
test.dt <- read.csv("00872872.csv")
den.dt <- density(test.dt$DNA_Index)
plot(den.dt)

var.est <- getVarianceViaDeriv (den.dt$x, den.dt$y)
var.est
mean.est <- getMeanViaDeriv(den.dt$x, den.dt$y)
mean.est




