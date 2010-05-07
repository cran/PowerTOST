# Author: dlabes
#------------------------------------------------------------------------------
# helper function to calculate std err from CV of lognormal data
.CV2se <- function(CV) return(sqrt(log(1.0 + CV^2)))
# reverse helper function
.se2CV <- function(se) return(sqrt(exp(se*se)-1))
#------------------------------------------------------------------------------
# 'raw' power function without any error checks,
# does not vectorize propperly!
# to be used in sampleN.TOST avoiding overhead of redundant calculations
# in case of multiplicative model:
# diffm=log(null ratio), theta1=log(lower BE limit), theta2=log(upper BE limit)
# in case of additive model:
# diffm=1-null ratio, theta1=lower BE limit-1, theta2=upper BE limit -1
.power.TOST <- function(alpha=0.05, theta1, theta2, diffm, se, n, df, bk=2)
{
  tval   <- qt(1 - alpha, df, lower.tail = TRUE)
  
  delta1 <- (diffm-theta1)/(se*sqrt(bk/n))
  delta2 <- (diffm-theta2)/(se*sqrt(bk/n))
  R      <- (delta1-delta2)*sqrt(df)/(2.*tval)
  
  # to avoid numerical errors in OwensQ implementation
  if (df[1]>10000) { 
    # Joulious formula (57) or (67), normal approximation
    p1 <- pnorm( (abs(delta1)-tval), lower.tail = TRUE)
    p2 <- pnorm( (abs(delta2)-tval), lower.tail = TRUE)
		
    return (p1 + p2 -1.)
  }
  # attempt to vectorize (it vectorizes properly if diffm is a vector
  # OR se is a vector OR n,df) TODO: check it!
  nel <- length(delta1)
  dl <- length(tval)
  p1 <-c(1:nel)	
  p2 <- p1
  for (i in seq_along(delta1)) {
	if (dl>1) {ddf <- df[i]; ttt <- tval[i]} 
	else {ddf <- df[1]; ttt <- tval[1]}
	p1[i] <- OwensQ(ddf,  ttt, delta1[i], 0, R[i])
	p2[i] <- OwensQ(ddf, -ttt, delta2[i], 0, R[i])
  }
  return( p2-p1 )
}
#------------------------------------------------------------------------------
# 'raw' power function without any error checks, 
# approximation based on non-central t
# this vectorices ok
.approx.power.TOST <- function(alpha=0.05, theta1, theta2, diffm, 
		                       se, n, df, bk=2)
{
  tval <- qt(1 - alpha, df, lower.tail = TRUE, log.p = FALSE)
  delta1 <- (diffm-theta1)/(se*sqrt(bk/n))
  delta2 <- (diffm-theta2)/(se*sqrt(bk/n))
	
  pow <- pt(-tval, df, ncp=delta2)-pt(tval, df, ncp=delta1)
  pow[pow<0] <- 0 # this is to avoid neg. power due to approx. (vector form)
  
  return(pow)
}
#------------------------------------------------------------------------------
# 'raw' power function without any error checks, 
# approximation based on central 'shifted' t distribution
# according to Chow, Liu "Design and Analysis of Bioavailability ..."
# Chapter 9.6 and implemented in PASS 2008
# where does this all come from?
.approx2.power.TOST <- function(alpha=0.05, theta1, theta2, diffm, 
                                se, n, df, bk=2)
{
	tval <- qt(1 - alpha, df, lower.tail = TRUE)
	delta1 <- (diffm-theta1)/(se*sqrt(bk/n))
	delta2 <- (diffm-theta2)/(se*sqrt(bk/n))
	
	pow <- pt(-delta2-tval,df) - pt(tval-delta1,df)
	pow[pow<0] <- 0 # this is to avoid neg. power due to approx. (vector form)
	
	return(pow)
}
#------------------------------------------------------------------------------
# Power of two-one-sided-t-tests using OwensQ or approx. using non-central t
# (this is a wrapper to .power.TOST(...) and .approx.power.TOST(...))
# In case of logscale=TRUE give diff, theata1 and theta2 as ratios
# f.i. ldiff=0.95, theta1=0.8, theta2=1.25
# In case of logscale=FALSE give diff, theata1 and theta2 as difference 
# to 1 f.i. diff=0.05 (5% difference), 
# theata1=-0.2, theta2=0.2 20% equiv. margins)
# CV is always the coefficient of variation but as ratio, not % 
# leave upper BE margin (ltheta2) empty and the function will use -lower
# in case of additive model or 1/lower if logscale=TRUE
power.TOST <- function(alpha=0.05, logscale=TRUE, theta1=0.8, theta2, 
		                   diff=0.95, CV, n, design="2x2", exact=TRUE)
{
  # check if design is implemented
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Err: unknown design!")
  
  # design characteristics
  ades <- .design.props(d.no)
  dfe  <- parse(text=ades$df[1],srcfile=NULL) #degrees of freedom as expression
  bk   <- ades$bk                             #design const.
  
  if (missing(CV)) stop("Err: CV must be given!")
  if (missing(n))  stop("Err: number of subjects must be given!")
  # handle log-transformation	
  if (logscale) {
    if (missing(theta2)) theta2 <- 1/theta1
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    ldiff   <- log(diff)
    se      <- sqrt(log(1.+CV^2))
  } else {
    if (missing(theta2)) theta2 <- -theta1
    ltheta1 <- theta1
    ltheta2 <- theta2
    ldiff   <- diff
    se      <- CV
  }
	
  df <- eval(dfe)
  if ( !exact )
    pow <- .approx.power.TOST(alpha, ltheta1, ltheta2, ldiff, se, n, df, bk)
  else
    pow <- .power.TOST(alpha, ltheta1, ltheta2, ldiff, se, n, df, bk)
	
  return( pow )
}

#------------------------------------------------------------------------------
# Approximate "expected" power according to Joulious book
# taking into account the uncertainty of an estimated se with 
# dfse degrees of freedom
# Only for log-transformed data.
# Raw function: see the call in exp.power.TOST()
.exppower.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, 
                            se, dfse, n, df, bk=2)
{
  tval <- qt(1 - alpha, df, lower.tail = TRUE)
  d1   <- sqrt(n*(diffm-ltheta1)^2/(bk*se^2))
  d2   <- sqrt(n*(diffm-ltheta2)^2/(bk*se^2))
  
  pow  <- pt(d1,dfse,tval) + pt(d2,dfse,tval) - 1
  
  return(pow)
  
}  
#------------------------------------------------------------------------------
# Joulious "expected" power
exppower.TOST <- function(alpha=0.05, theta1=0.8, theta2, diff=0.95, 
                           CV, dfCV, n, design="2x2")
{
  # check if design is implemented
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Err: unknown design!")
  
  # design characteristics
  ades <- .design.props(d.no)
  dfe  <- parse(text=ades$df[1],srcfile=NULL) #degrees of freedom as expression
  bk   <- ades$bk                             #design const.
  
  if (missing(CV) | missing(dfCV)) stop("Err: CV and/or df must be given!")
  if (missing(n))  stop("Err: number of subjects must be given!")

  if (missing(theta2)) theta2 <- 1/theta1
  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  ldiff   <- log(diff)
  se      <- sqrt(log(1.+CV^2))
  df <- eval(dfe)
  pow <- .exppower.TOST(alpha, ltheta1, ltheta2, ldiff, se, dfCV, n, df, bk)

  return( pow )
}