#-- The functions of normal-, t-distributions and integrate() ------------------
#require(stats) #this is usually not necessary within a standard installation
#-------------------------------------------------------------------------------
# Owen's Q-function 
# a, b must be a scalar numeric
# nu, t and delta also, no vectors allowed
OwensQ <- function (nu, t, delta, a, b)
{
	if (length(nu)>1 | length(t)>1 | length(delta)>1 | length(a)>1 | 
      length(b)>1) stop("Input must be scalars!")
  if (nu<1) stop("nu must be >=1!")
  # Observation: for 'really' large df (nu>2000) and large delta/b the  
	# density function is zero over nearly all its range! Q than returned 
	# sometimes falsly as =0! See documentation ?integrate for that.
	# example: OwensQ(3000,1.64,-10,0,300) = 1
	#          OwensQ(3000,1.64,-10,0,350) = 5.614476e-12 !
	# Idea: adapt upper and/or lower integration limit to account for that
	low <- a; up <- b
  # May 2011: shrink  interval depending on delta*b 
  if (b>1E6) b <- Inf
  if (is.finite(b)){ # in case of alpha=0.5 b is infinite
  	if (nu >= 1000 || abs(delta*b) > 30 || b>50){
      # try to shorten the range via interval halving: Jul 2012
      # upper integration limit
      ab <- b-a
      x  <- b
      dens <- .Q.integrand(x, nu, t, delta)
      #cat("Upper search\n")
      #cat("x=",x,"dens=",dens,"\n")
      while (dens==0){
        ab   <- ab*0.5
        x    <- x - ab
        dens <- .Q.integrand(x, nu, t, delta)
        #cat("x=",x,"dens=",dens,"\n")
        if (ab < 1E-10) break
      }
      if (dens>0) up <- min(x + ab, b) else return(0)
      #cat("up=",up,"\n")
      # lower limit
      ab   <- up - a
      x    <- a
      dens <- .Q.integrand(x, nu, t, delta)
      #cat("Lower search\n")
      #cat("x=",x,"dens=",dens,"\n")
      while (dens==0){
        ab   <- ab*0.5 
        x    <- x + ab
        dens <- .Q.integrand(x, nu, t, delta)
        #cat("x=",x,"dens=",dens,"\n")
        if (ab<1E-10) break
      }
      if (dens>0)  low <- max(x - ab, a) else return(0)
      #cat("low=",low,"\n")
    }  
  }
	# result of integrate() is a list, see ?integrate
	# .Machine$double.eps^.5 = 1.490116e-08 on my machine
	# MBESS uses .Machine$double.eps^0.25 = 0.0001220703 for both tolerances
	# seems it makes no difference
	Qintegral <- integrate(.Q.integrand, lower = low, upper = up, 
			         nu=nu, t=t, delta = delta, subdivisions = 10000, 
			         #rel.tol = .Machine$double.eps^0.5, 
               #abs.tol = .Machine$double.eps^0.5,
			         rel.tol = 1.e-10, abs.tol=1.e-12, stop.on.error = TRUE)
	# error handling? How?
	return(Qintegral[[1]])  
}
#-------------------------------------------------------------------------------
# Integrand of the definit integral in Owen's Q. Used in the call of integrate()
# Not useful alone, I think ? Leading . hides this function 
# function must give a vectorized answer in respect to x
.Q.integrand <- function(x, nu, t, delta)
{ #version without for - loop, it works without
	lnQconst <- -((nu/2.0)-1.0)*log(2.0) - lgamma(nu/2.)

# what if x<0? Should here not possible, but ...
# simple x^(nu-1) doesnt work for high nu because  = Inf 
# and then exp( -0.5*x^2 + lnQconst )*x^(nu-1) -> NaN
# (nu-1)*log(abs(x)) is NaN if nu=1, x=0! 0*(-Inf) -> NaN
  
  dens <- x  # assures that dens=0 if x=0
  dens[x!=0] <- sign(x)^(nu-1) *
      pnorm( t*x/sqrt(nu) - delta, mean = 0, sd = 1, log.p = FALSE) * 
      exp( (nu-1)*log(abs(x)) - 0.5*x^2 + lnQconst )
      
	dens
}

# Test cases:
# Craig Zupke's observations:
# power.TOST(0.410,FALSE,-5.97,5.97,8.5448,1,14,"parallel",TRUE) #!old call
# power.TOST(0.410,FALSE,-5.97,5.97,8.5448,1,14,"parallel","exact")
# gave an error; high b/delta
# should give: 2.335633e-07

# Jul 2012: Helmuts observation
# n=4, CV=1E-5(=se) gives power=1 (delta1=24303.3, delta2=-38811.23, R=b=15283.88
#      CV=1E-6 gives power=0      (      243033          -388112.3   R  152838.8
#      CV=0    gives power=1             Inf              -Inf       Inf
# tval=2.919986
# for CV=1e-6: erroneous in versions pre 0.9-9. 2. call gave =0
# OwensQ(nu=2, t= 2.919986, delta= 243033,   0, 152838.8) ==0
# OwensQ(nu=2, t=-2.919986, delta=-388112.3, 0, 152838.8) ==1
# for CV=0
# OwensQ(nu=2, t=2.919986, delta=Inf, 0, Inf)  ==0
# OwensQ(nu=2, t=-2.919986, delta=-Inf, 0,Inf) ==1
# 
# Helmuts cases (ver 0.9-9) Jul 2012
# sampleN.TOST(theta0=1, CV=0.02, design="2x2", print=TRUE) # Ok
# #next gave an error due to 0*-Inf in .Q.integrand()
# sampleN.TOST(theta0=1, CV=0.01, design="2x2", print=TRUE) 

