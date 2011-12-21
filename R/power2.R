###############################################################################
# Function(s) for power calculation for unbalanced (sequence) groups
# 
# Author: dlabes
###############################################################################

power2.TOST <- function(alpha=0.05, logscale=TRUE, theta1, theta2, theta0,
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
  
  nc <- n # for use in call of the raw functions since n is overwritten
  bk <- ades$bkni
  steps <- ades$steps
  if (is.na(bk)) stop("Not able to handle unbalanced ",design," design!")
  if (design=="parallel") {
    dfe   <- parse(text="n-2",srcfile=NULL)# for total
    steps <- 2
  }
  # check the correct length due to design
  # must be length==number of (sequence) groups, that is coded in steps
  # except for the parallel group design
  if (length(n)!= steps){
     stop("Length of n vector must be ",steps,"!")
  }
  nc <- sum(1/n)  # 1/nc is used in call of the raw functions
                  # with bkni this results in bkni/(1/nc)=bkni*sum(1/n)
  n  <- sum(n)    # for use in df expression
  
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
  if (any(df<1)) stop("n's too small. Degrees of freedom <1!")
  
  pow <- .calc.power(alpha, ltheta1, ltheta2, ldiff, se, 1/nc, df, bk, method)
  
  return(pow)
}
