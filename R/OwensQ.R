#-- The functions of normal-, t-distributions and integrate() ------------------
require(stats)
#-------------------------------------------------------------------------------
# Owen's Q-function 
OwensQ <- function (nu, t, delta, a, b)
{
	# Observation: for 'really' large df (nu>2000) and large delta/b the  
	# density function is zero over nearly all its range. Q is than returned 
	# sometimes as =0! See documentation ?integrate for that.
	# example: OwensQ(3000,1.64,-10,0,300) = 1
	#          OwensQ(3000,1.64,-10,0,350) = 5.614476e-12 !
	# Idea: adapt upper and/or lower to account for that
	low <- a; up<- b
	if (nu>=1000){
		# try to shorten the range over that is integrated
		x <- seq(a, b, by=(b-a)/400)
		dens <- .Q.integrand(x, nu, t, delta)
		r <- data.frame(x=x, dens=dens)
		r <- r[r$dens>0,]
		if (length(r)>0) {# if any >0
			up <- max(r$x)+(b-a)/600  # only the upper x+step is used
			rm(r,x,dens)
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
	dens <- pnorm( t*x/sqrt(nu) - delta, mean = 0, sd = 1, log = FALSE) * 
			exp( (nu-1)*log(x) - 0.5*x^2 + lnQconst )
	
	return(dens)
}