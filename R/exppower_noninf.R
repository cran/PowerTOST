#------------------------------------------------------------------------------
# Author: dlabes
#------------------------------------------------------------------------------
# Approximate, "expected" power according to Julious book
# taking into account the uncertainty of an estimated se with 
# dfse degrees of freedom
# Only for log-transformed data.
# Raw function: see the call in exppower.noninf()
.exppower.noninf <- function(alpha=0.05, lmargin, diffm, sedm, dfse, df)
{
  tval <- qt(1 - alpha, df, lower.tail = TRUE)
  tau  <- sqrt((diffm-lmargin)^2/sedm^2)
  # in case of diffm=lmargin and se=0
  # tau has the value NaN
  tau[is.nan(tau)] <- 0
  
  pow  <- pt(tau,dfse,tval)
  
  return(pow)
  
}  
#------------------------------------------------------------------------------
# Julious "expected" power
exppower.noninf <- function(alpha=0.025, logscale=TRUE, theta0, margin, 
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
  bk   <- ades$bk
  
  if (missing(CV) | missing(dfCV)) stop("CV and df must be given!", call.=FALSE)
  if (missing(n)) stop("Number of subjects must be given!", call.=FALSE)
  
  if (logscale){
    if (missing(theta0)) theta0 <- 0.95
    if (missing(margin)) margin <- 0.8
    lmargin <- log(margin)
    ldiff   <- log(theta0)
    se      <- CV2se(CV)
  } else {
    if (missing(theta0)) theta0 <- -0.05
    if (missing(margin)) margin <- -0.2
    lmargin <- margin
    ldiff   <- theta0
    se      <- CV
  }
  
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
  pow <- .exppower.noninf(alpha, lmargin, ldiff, se*se.fac, dfCV, df)
  
  return( pow )
}



