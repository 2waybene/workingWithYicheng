##====================================================
##	File:   simulate_DI.R
##	Author: Jianying Li
##	Comment: Learning mixture model for Aneuploidy
##		   study for OSCC
##====================================================
root <- "x:/"
root <- "/Users/li11/"
source( paste (root, "myGit/workingWithYicheng/utils/", "mixtureModelFunctions.R", sep = ""))
##=====End of settings
##=====Sample data
setwd( paste (root, "myGit/workingWithYicheng/sampleData", sep=""))


#f <- "fmd_DT_raw.txt"


##==================================================
##	Getting coefficient of variation 
##	from a real data
##==================================================
f <- "oscc-olk1_parsed.txt"
f_IN <- paste (root, "myGit/workingWithYicheng/sampleData/",  f , sep="")
dt <- read.table(f_IN, header= F, sep = "\t")
str(dt)
temp.cv <- sd(dt$V1)/mean(dt$V1)
temp.cv

##=====End

##======================================================================
##  Based on our prior knowledge, we can focus on 
##  three clusters of cell populations:
##  Cluster one, with mean DI = 1.001
##  Cluster two, with mean DI = 2.001
##  Cluster three, with minimum DI > 2.300, assigning 2.300 as mean
##  Now, let's assume that C.V. is same 
##======================================================================
temp.cv

mean1 <- c(1.001, 2.002, 2.300)
mean2 <- c(1.001, 2.002, 3.600)
mean3 <- c(0.936, 1.692, 3.600)

weight1 <- c(0.893, 0.092, 0.05)
weight2 <- c(0.425, 0.425, 0.05)

plot_Mixed_3_families <- function (mean, sigma , weight)
{

	Delta <- 0.01
	x <- seq(0,7, by=Delta)
	y1 <- weight[1]*P(x, mean[1], sigma[1])
	y2 <- weight[2]*P(x, mean[2], sigma[2])
	y3 <- weight[3]*P(x, mean[3], sigma[3])
	y <- y1 + y2 + y3
	
	par(mfrow=c(1,1))
	plot(x,y, type="l", lwd=3,
     		main="Simulated D.I. values: three-category groups",
	     xlab="D.I. value", ylab="Probability Density")
	abline(h=0, lty="dotted")
	lines(x,y1, col="red")
	lines(x,y2, col="green")
	lines(x,y3, col="blue")

	lgd = c("Mixture", paste("Normal, mean = ", mean[1], sep=""), 
			paste ("Mitotic,  mean = ",  mean[2], sep= ""), 
			paste ("Aneuploid, mean = ", mean[3], sep= ""))
	legend ("topright", lgd, text.col = c("black", "red", "green", "blue"))
}

sigma1 <- c(0.19, 0.25, 0.5) 
sigma2 <- mean1*temp.cv

##	get a few plots..

plot_Mixed_3_families (mean1, sigma1, weight1)
plot_Mixed_3_families (mean1, sigma2, weight1)
plot_Mixed_3_families (mean1, sigma1, weight2)
plot_Mixed_3_families (mean3, sigma1, weight2)

##=======================================================
##	Estimating the mean and standard deviation from
##	derivatives
##=======================================================

plot_derivatives <- function (mean, sigma , weight)
{

#mean   <-  mean1
#sigma  <-  sigma1
#weight <-  weight1

	Delta <- 0.01
	x <- seq(0,7, by=Delta)
	y1 <- weight[1]*P(x, mean[1], sigma[1])
	y2 <- weight[2]*P(x, mean[2], sigma[2])
	y3 <- weight[3]*P(x, mean[3], sigma[3])
	y <- y1 + y2 + y3

	derivative1 <- Deriv1(x,y)
	derivative2 <- Deriv2(x,y)

	par(mfrow=c(3,1))
	plot(x,y, type="l", lwd=3,
		main="Simulated D.I. values: three-category groups",
	  	xlab="D.I. value", ylab="Probability Density")
		abline(h=0, lty="dotted")
	lines(x,y1, col="red")
	lines(x,y2, col="green")
	lines(x,y3, col="blue")

	##	2nd figure
	##### 1st Derivative

	plot(derivative1$x,derivative1$y, type="l",
		main="1st Derivative", xlab="x", ylab="y'")
		abline(h=0, lty="dotted")

	peaks.Deriv1   <- peaks( derivative1$y, span=3)
	valleys.Deriv1 <- peaks(-derivative1$y, span=3)

	points( derivative1$x[peaks.Deriv1], derivative1$y[peaks.Deriv1],
      	  pch=19, col="red")
	points( derivative1$x[valleys.Deriv1], derivative1$y[valleys.Deriv1],
      	  pch=19, col="blue")

	# Approximate location of peak and valley
	derivative1$x[peaks.Deriv1]
	derivative1$x[valleys.Deriv1]

	s.approx <- (derivative1$x[valleys.Deriv1] - derivative1$x[peaks.Deriv1])/2
	
	#print(s.approx)

	# Approximate standard deviation
	print ("Here are the standard deviation")
	sd.temp <- (derivative1$x[valleys.Deriv1] -  derivative1$x[peaks.Deriv1])/2
	print (sd.temp)
	##### 2nd Derivative

	plot(derivative2$x,derivative2$y, type="l",
		main="2nd Derivative", xlab="x", ylab="y''")
		abline(h=0, lty="dotted")

	peaks.Deriv2   <- peaks( derivative2$y, span=3)
	valleys.Deriv2 <- peaks(-derivative2$y, span=3)

	points( derivative2$x[peaks.Deriv2], derivative2$y[peaks.Deriv2],
      	  pch=19, col="red")
	points( derivative2$x[valleys.Deriv2], derivative2$y[valleys.Deriv2],
      	  pch=19, col="blue")

	# Approximate location of mean of normal distribution:
	derivative2$x[valleys.Deriv2]
	derivative2$y[valleys.Deriv2]

	print ("Here are the means")
	mu.approx <-  derivative2$x[valleys.Deriv2]
	print(mu.approx)
}

##==========================================================
plot_derivatives (mean1, sigma1, weight1)
plot_derivatives (mean1, sigma2, weight1)
plot_derivatives (mean1, sigma1, weight2)
plot_derivatives (mean3, sigma1, weight2)



##===========================================================
##	old codes
##	will revisit later: JYL
##	11/06/2015
##====================
	# Peaks
derivative2$x[peaks.Deriv2]
derivative2$y[peaks.Deriv2]


# Valleys
derivative2$x[valleys.Deriv2]
derivative2$y[valleys.Deriv2]



# Attempt non-linear curve fit to extract the parameters
set.seed(17)
y <- y + rnorm(length(y), 1E-5, 1E-4)

fit.pike <- nls(y ~
                (a/b)*exp(-(x-c)^2/(2*b^2)) +
                (d/e)*exp(-(x-f)^2/(2*e^2)),
 
                start=list(a=(1/sqrt(2*pi)) / mean[1], b=mean[1], c=sigma[1],
                           d=(1/sqrt(2*pi)) / mean[2], e=mean[2], f=sigma[2]),
                control=nls.control(tol=1E-5, minFactor=1/1024),
        	trace=TRUE)

# Means (mu values)
coef(fit.pike)[3*1:4]

# Standard deviations (s values)


 
##========================================================
## Summary data from 29 samples
##========================================================


f <- "summary_dt.txt";
f_IN <- paste (file.dir, f , sep="")
f_IN <- paste (macFileDir, f , sep="")
f_IN <- paste (winFileDir, f , sep="")
##=====End of settings

dt <- read.table(f_IN, header= T, sep = "\t")

str(dt)

summary(dt$Mean3)



