#-- The functions of normal-, t-distributions and integrate() ------------------
require(stats)
#-------------------------------------------------------------------------------
# Owen's Q-function 
# a, b must be a scalar numeric
# nu, t and delta also, no vectors allowed
OwensQ <- function (nu, t, delta, a, b)
{
	# Observation: for 'really' large df (nu>2000) and large delta/b the  
	# density function is zero over nearly all its range! Q is than returned 
	# sometimes falsly as =0! See documentation ?integrate for that.
	# example: OwensQ(3000,1.64,-10,0,300) = 1
	#          OwensQ(3000,1.64,-10,0,350) = 5.614476e-12 !
	# Idea: adapt upper and/or lower integration limit to account for that
	low <- a; up <- b
	if (nu>=1000){
		# try to shorten the integration range
    h <- (b-a)/499
		x <- seq(a, b, by=h);x[1] <- a; x[500] <- b
		dens <- .Q.integrand(x, nu, t, delta)
		x <- x[dens > 0] # .Machine$double.eps approx. 2.22e-16
		n <- length(x)
    if (n > 0) {# if any >0
			up  <- min(x[n]+h, b)  # upper: xoff+step if < b
			low <- max(x[1]-h, a)  # lower  xon-step if > a
		} else {
      # all == 0
      return(0.0)
    }
	}
	# result of integrate() is a list, see ?integrate
	# .Machine$double.eps^.5 = 1.490116e-08 on my machine
	# MBESS uses .Machine$double.eps^0.25 = 0.0001220703
	# seems it makes no difference
	Qintegral <- integrate(.Q.integrand, lower = low, upper = up, 
			         nu=nu, t=t, delta = delta, subdivisions = 10000, 
			         #rel.tol = .Machine$double.eps^0.5, 
			         #abs.tol = .Machine$double.eps^0.5,
			         rel.tol = 1.e-10, abs.tol=1.e-10, 
			         stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)
	Q <- Qintegral[[1]]
	
	# error handling? How?
	return(Q)  
}
#-------------------------------------------------------------------------------
# Integrand of the definit integral in Owen's Q. Used in the call of integrate()
# Not useful alone, I think ? Leading . hides this function    
.Q.integrand <- function(x, nu, t, delta)
{ #version without for - loop, it works without
	lnQconst <- -((nu/2.0)-1.0)*log(2.0) - lgamma(nu/2.)
# original code in first release:
#  dens <- pnorm( t*x/sqrt(nu) - delta, mean = 0, sd = 1, log = FALSE) * 
#          exp( (nu-1)*log(x) - 0.5*x^2 + lnQconst )
# what if x<0? Should here not possible, but ...
# simple x^(nu-1) doesnt work for high nu because  = inf 
# and then exp( -0.5*x^2 + lnQconst )*x^(nu-1) -> NaN
  dens <- sign(x)^(nu-1) *
          pnorm( t*x/sqrt(nu) - delta, mean = 0, sd = 1, log = FALSE) * 
          exp( (nu-1)*log(abs(x)) - 0.5*x^2 + lnQconst )
      
	return(dens)
}