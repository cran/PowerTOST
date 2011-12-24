# Author: dlabes
#------------------------------------------------------------------------------
# helper function to calculate std err from CV of lognormal data
CV2se <- function(CV) return(sqrt(log(1.0 + CV^2)))
# reverse helper function
se2CV <- function(se) return(sqrt(exp(se*se)-1))
#------------------------------------------------------------------------------
# helper function to allow partial match of the method
.powerMethod <- function(method){
  meth <- tolower(method[1])
  if (method=="") method <- "exact"
  methods <- c("exact","owenq","noncentral","nct","shifted","central")
  #                         ^ = match at start
  meth <- methods[grep(paste("^",meth,sep=""), methods)]
  if (length(meth)==0) meth <- tolower(method)
  if (length(meth)>1)  meth <- "nct" # happens if only "n" is given
  
  return(meth)
}
#------------------------------------------------------------------------------
# 'raw' power function without any error checks,
# does not vectorize propperly!
# to be used in sampleN.TOST avoiding overhead of redundant calculations
# in case of multiplicative model:
# diffm=log(null ratio), theta1=log(lower BE limit), theta2=log(upper BE limit)
# in case of additive model:
# diffm=1-null ratio, theta1=lower BE limit-1, theta2=upper BE limit -1
.power.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, se, n, df, bk=2)
{
  tval   <- qt(1 - alpha, df, lower.tail = TRUE)
  # if alpha>0.5 (very unusual) then R is negative 
  # in the application of OwensQ the upper integration limit 
  # is lower then the lower integration limit!
  # SAS OwenQ gives missings if b or a are negative!
  
  delta1 <- (diffm-ltheta1)/(se*sqrt(bk/n))
  delta2 <- (diffm-ltheta2)/(se*sqrt(bk/n))
  # R is infinite in case of alpha=0.5
  R      <- (delta1-delta2)*sqrt(df)/(2.*abs(tval))
  
  # to avoid numerical errors in OwensQ implementation
  if (df[1]>50000) { 
    # Joulious formula (57) or (67), normal approximation
    p1 <- pnorm( (abs(delta1)-tval), lower.tail = TRUE)
    p2 <- pnorm( (abs(delta2)-tval), lower.tail = TRUE)
		
    return(p1 + p2 -1.)
  }
  if (df[1]>10000 & df[1]<=50000) {
    # approximation via non-central t-distribution
    return(.approx.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk))
  }
  
  # attempt to vectorize (it vectorizes properly if diffm is a vector
  # OR se OR n,df are vectors) 
  nel <- length(delta1)
  dl <- length(tval)
  p1 <- c(1:nel)	
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
# 'raw' approximate power function without any error checks, 
# approximation based on non-central t
# this vectorizes ok
.approx.power.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, 
		                           se, n, df, bk=2)
{
  tval <- qt(1 - alpha, df, lower.tail = TRUE, log.p = FALSE)
  delta1 <- (diffm-ltheta1)/(se*sqrt(bk/n))
  delta2 <- (diffm-ltheta2)/(se*sqrt(bk/n))
	
  pow <- pt(-tval, df, ncp=delta2)-pt(tval, df, ncp=delta1)
  pow[pow<0] <- 0 # this is to avoid neg. power due to approx. (vector form)
  
  return(pow)
}
#------------------------------------------------------------------------------
# 'raw' power function without any error checks, 
# approximation based on central 'shifted' central t distribution
# according to Chow, Liu "Design and Analysis of Bioavailability ..."
# Chapter 9.6 and implemented in PASS 2008
# where does this all come from?
.approx2.power.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, 
                                se, n, df, bk=2)
{
	tval <- qt(1 - alpha, df, lower.tail = TRUE)
  delta1 <- (diffm-ltheta1)/(se*sqrt(bk/n))
	delta2 <- (diffm-ltheta2)/(se*sqrt(bk/n))
	
	pow <- pt(-delta2-tval,df) - pt(tval-delta1,df)
	pow[pow<0] <- 0 # this is to avoid neg. power due to approx. (vector form)
	
	return(pow)
}
# function merging the various power calculations
.calc.power <- function(alpha=0.05, ltheta1, ltheta2, diffm, se, n, df, bk, 
                        method="exact")
{
  pow <- switch(
      method,
      exact=.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk),
      owenq=.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk),
      nct=  .approx.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk),
      noncentral=.approx.power.TOST(alpha, ltheta1, ltheta2, diffm, 
                                    se, n, df, bk),
      shifted=.approx2.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk),
      central=.approx2.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk),
      stop("Method ", method, " unknown!", .call=TRUE)
  ) 
  return(pow)
}
#------------------------------------------------------------------------------
# Power of two-one-sided-t-tests using OwensQ or approx. using non-central t
# (this is a wrapper to .power.TOST(...) and .approx.power.TOST(...))
# In case of logscale=TRUE give theta0, theata1 and theta2 as ratios
# f.i. theta0=0.95, theta1=0.8, theta2=1.25
# In case of logscale=FALSE give theta0, theata1 and theta2 as difference 
# to 1 f.i. theta0=0.05 (5% difference), 
# theata1=-0.2, theta2=0.2 20% equiv. margins)
# CV is always the coefficient of variation but as ratio, not % 
# leave upper BE margin (ltheta2) empty and the function will use -lower
# in case of additive model or 1/lower if logscale=TRUE
power.TOST <- function(alpha=0.05, logscale=TRUE, theta1, theta2, theta0,
                       CV, n, design="2x2", method="exact", robust=FALSE)
{
  # check if design is implemented
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design, " unknown!", call.=FALSE)
  
  # design characteristics
  ades <- .design.props(d.no)
  #degrees of freedom as expression
  if (robust){
    dfe  <- parse(text=ades$df[1],srcfile=NULL) 
  } else {
    dfe  <- parse(text=ades$df2[1],srcfile=NULL)
  }
  if (missing(CV)) stop("CV must be given!")
  if (missing(n))  stop("Number of subjects n must be given!")
  
  #design const.
  bk    <- ades$bk
  steps <- ades$steps
  
  # regularize the method giving
  method <- .powerMethod(method)
  
  # handle log-transformation	
  if (logscale) {
    if (missing(theta0)) theta0 <- 0.95
    if (missing(theta1)) theta1 <- 0.8
    if (missing(theta2)) theta2 <- 1/theta1
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    ldiff   <- log(theta0)
    se      <- CV2se(CV)
  } else { # untransformed
    if (missing(theta1)) theta1 <- -0.2
    if (missing(theta0)) theta0 <- 0.05
    if (missing(theta2)) theta2 <- -theta1
    ltheta1 <- theta1
    ltheta2 <- theta2
    ldiff   <- theta0
    se      <- CV
  }
  
  df <- eval(dfe)
  if (any(df<1)) stop("n too small. Degrees of freedom <1!")
  pow <- .calc.power(alpha, ltheta1, ltheta2, ldiff, se, n, df, bk)
  return( pow )
}
