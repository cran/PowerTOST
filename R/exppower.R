#------------------------------------------------------------------------------
# Author: dlabes
#------------------------------------------------------------------------------
# Approximate "expected" power according to Julious book
# taking into account the uncertainty of an estimated se with 
# dfse degrees of freedom
# Only for log-transformed data.
# Raw function: see the call in exppower.TOST()
.exppower.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, sem, dfse, df)
{
  tval <- qt(1 - alpha, df, lower.tail = TRUE)
  d1   <- sqrt((diffm-ltheta1)^2/sem^2)
  d2   <- sqrt((diffm-ltheta2)^2/sem^2)
  # in case of diffm=ltheta1 or =ltheta2 and se=0
  # d1 or d2 have then value NaN (0/0)
  d1[is.nan(d1)] <- 0
  d2[is.nan(d2)] <- 0
  
  pow  <- pt(d1,dfse,tval) + pt(d2,dfse,tval) - 1
  
  return(pow)
  
}  
#------------------------------------------------------------------------------
# Julious "expected" power
exppower.TOST <- function(alpha=0.05, logscale=TRUE, theta0, theta1, theta2, 
                          CV, dfCV, n, design="2x2", robust=FALSE)
{
  # check if design is implemented
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design, " unknown!", call.=FALSE)
  
  # design characteristics
  ades <- .design.props(d.no)
  #df as expression
  dfe  <- .design.df(ades, robust=robust)
  #design const.
  #bk   <- ades$bk # we use always bkni
  
  if (missing(CV) | missing(dfCV)) stop("CV and df must be given!", call.=FALSE)
  if (missing(n)) stop("Number of subjects must be given!", call.=FALSE)
  
  if (logscale) {
    if (missing(theta0)) theta0 <- 0.95
    if (missing(theta1)) theta1 <- 0.8 
    if (missing(theta2)) theta2 <- 1/theta1
    
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    ldiff   <- log(theta0)
    se      <- CV2se(CV)
    
  } else {
    if (missing(theta0)) theta0 <- 0.05
    if (missing(theta1)) theta1 <- -0.2 
    if (missing(theta2)) theta2 <- -theta1
    ltheta1 <- theta1
    ltheta2 <- theta2
    ldiff   <- theta0
    se      <- CV
  }
  # we use always bkni
  #bk   <- ades$bk
  if (length(n) == 1) {
    # total n given    
    # for unbalanced designs we divide the ns by ourself
    # to have only small imbalance (function nvec() from Helper_dp.R)
    n <- nvec(n=n, grps=ades$steps)
    if (n[1]!=n[length(n)]){
      message("Unbalanced design. n(i)=", paste(n, collapse="/"), " assumed.")
    } 
  } else {
    if (length(n) != ades$steps) {
      stop("Length of n vector must be ", ades$steps, "!")
    }
  }
  
  nc <- sum(1/n)
  n <- sum(n)
  se.fac <- sqrt(ades$bkni * nc)
  
  
  df      <- eval(dfe)
  if (any(df<1)) stop("n too low. Degrees of freedom <1!", call.=FALSE)
  
  pow <- .exppower.TOST(alpha=alpha, ltheta1=ltheta1, ltheta2=ltheta2, 
                        diffm=ldiff, sem=se*se.fac, dfse=dfCV, df=df)
  
  return( pow )
}



