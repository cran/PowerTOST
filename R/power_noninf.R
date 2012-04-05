###############################################################################
# Power calculations based on non-inferiority t-test
# 
# Author: dlabes
###############################################################################

# --------------------------------------------------------------------------
# internal functions:
# power function (working horse)
.power.noninf <- function(alpha=0.025, lmargin, diffm, se, n, df, bk=2)
{
  tval <- qt(1-alpha,df)
  tau  <- abs( (diffm-lmargin)*sqrt(n)/sqrt(bk*se^2) )
  return(1 - pt(tval, df, tau))
}

# --------------------------------------------------------------------------
# Power function for non-inferiority t-test (OOST)
power.noninf <- function(alpha=0.025,  logscale=TRUE, margin, theta0, CV, n, 
                         design="2x2", robust=FALSE)
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
  bk <- ades$bk
  
  if (missing(CV)) stop("CV must be given!")
  if (missing(n))  stop("Number of subjects n must be given!")
  
  df   <- eval(dfe)    
  if (any(df<1)) stop("n too small. Degrees of freedom <1!")
  
  # handle log-transformation
  if (logscale) {
    if (missing(margin)) margin <- 0.8
    if (missing(theta0)) theta0 <- 0.95
    lmargin <- log(margin)
    diffm   <- log(theta0)
    se      <- CV2se(CV)
  } else {
    if (missing(margin)) margin <- -0.2
    if (missing(theta0)) theta0 <- -0.05
    lmargin <- margin
    diffm   <- theta0
    se      <- CV
  }
  return(.power.noninf(alpha, lmargin, diffm, se, n, df, bk))
}
