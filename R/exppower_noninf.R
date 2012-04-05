# Author: dlabes
#------------------------------------------------------------------------------
# Approximate, "expected" power according to Julious book
# taking into account the uncertainty of an estimated se with 
# dfse degrees of freedom
# Only for log-transformed data.
# Raw function: see the call in exppower.TOST()
.exppower.noninf <- function(alpha=0.05, lmargin, diffm, 
                             se, dfse, n, df, bk=2)
{
  tval <- qt(1 - alpha, df, lower.tail = TRUE)
  tau  <- sqrt(n*(diffm-lmargin)^2/(bk*se^2))
  
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
  if (robust){ 
    # 'robust' evaluation df's
    dfe  <- parse(text=ades$df2[1], srcfile=NULL) 
  } else {
    dfe  <- parse(text=ades$df[1], srcfile=NULL) 
  }
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
  df      <- eval(dfe)
  if (any(df<1)) stop("n too low. Degrees of freedom <1!", call.=FALSE)
  pow <- .exppower.noninf(alpha, lmargin, ldiff, se, dfCV, n, df, bk)
  
  return( pow )
}



