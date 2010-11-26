# Author: dlabes
#------------------------------------------------------------------------------
# Approximate "expected" power according to Julious book
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
# Julious "expected" power
exppower.TOST <- function(alpha=0.05, theta1=0.8, theta2, theta0=0.95,
                          CV, dfCV, n, design="2x2")
{
  # check if design is implemented
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design, " unknown!",call.=FALSE)
  
  # design characteristics
  ades <- .design.props(d.no)
  dfe  <- parse(text=ades$df[1],srcfile=NULL) #degrees of freedom as expression
  bk   <- ades$bk                             #design const.
  
  if (missing(CV) | missing(dfCV)) stop("CV and df must be given!",call.=FALSE)
  if (missing(n)) stop("Number of subjects must be given!",call.=FALSE)
  
  if (missing(theta2)) theta2 <- 1/theta1
  
  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  ldiff   <- log(theta0)
  se      <- CV2se(CV)
  df      <- eval(dfe)
  pow <- .exppower.TOST(alpha, ltheta1, ltheta2, ldiff, se, dfCV, n, df, bk)
  
  return( pow )
}



