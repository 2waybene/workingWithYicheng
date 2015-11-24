##======================================================
##	peak-functions.R
##	Here are the peak function collections
##	For EdTAR projects
##======================================================


# Here P is the same as dnorm (probability density function for normal
# distribution), but other functions could be tried here.
P <- function(x, mean, sd)
{
  variance <- sd^2
  exp(-(x-mean)^2/(2*variance)) / sqrt(2*pi*variance)
}


# integrate(P, -1, 1, mean=0, sd=1) is same as integrate(dnorm, -1, 1, mean=0, sd=1)
# integrate(P,-5,5, mean=0, sd=1)   # should be close to 1.0

# Find "peaks" in array.
# R equivalent of Splus peaks() function.
# http://finzi.psych.upenn.edu/R/Rhelp02a/archive/33097.html
# (see efg's posting to R-Help on 8 Feb 2007 about problem with ties.)
#
# peaks(c(1,4,4,1,6,1,5,1,1),3)
# [1] FALSE FALSE  TRUE FALSE  TRUE FALSE  TRUE

peaks <- function(series,span=3)
{
  z <- embed(series, span)
  s <- span%/%2
  v <- max.col(z, "first") == 1 + s   # take first if a tie
  result <- c(rep(FALSE,s),v)
  result <- result[1:(length(result)-s)]
  result
}


# First derivative.  Adjust x values to be center of interval.
# Spacing of x-points need not be uniform
Deriv1 <- function(x,y)
{
  y.prime <- diff(y) / diff(x)
  x.prime <- x[-length(x)] + diff(x)/2
  list(x = x.prime,
       y = y.prime)
}

# "Centered" 2nd-derivative. Spacing of x-points assumed to be uniform.
Deriv2 <- function(x,y)
{
  h <- x[2] - x[1]
  Range <- 2:(length(x)-1)  # Drop first and last points
  list(x = x[Range],
       y = (y[Range+1] - 2*y[Range] + y[Range-1]) / h^2)
}


getVarianceViaDeriv <- function (x, y)
{
	##### 1st Derivative
	# According to Abramowitz & Stegun, inflection points are at +/- sigma.
	
	derivative1 <- Deriv1(x,y)
	peaks.Deriv1   <- peaks( derivative1$y, span=3)
	valleys.Deriv1 <- peaks(-derivative1$y, span=3)

	s.approx <- (derivative1$x[valleys.Deriv1] - derivative1$x[peaks.Deriv1])/2
	#	This is standard deviation!!
	return(s.approx)

}


getMeanViaDeriv <- function (x, y)
{
	#### 2nd Derivative
	derivative2 <- Deriv2(x,y)

	peaks.Deriv2   <- peaks( derivative2$y, span=3)
	valleys.Deriv2 <- peaks(-derivative2$y, span=3)

	mu.approx <-  derivative2$x[valleys.Deriv2]
	return(mu.approx)
}




