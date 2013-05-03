#------------------------------------------------------------------------------
# function to calculate the CV from given 90% confidence interval
# 
# Author: dlabes
#------------------------------------------------------------------------------

CVfromCI <- function(point, lower, upper, n, design="2x2", alpha=0.05, 
                     robust=FALSE)
{
  if (missing(lower) | missing(upper)) {
    stop("Lower and upper CL must be given!", call.=FALSE) 
  }
  if (missing(n)) stop("Sample size n must be given!",call.=FALSE)
  
  if (missing(point)) point <- sqrt(lower*upper) 

  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design," unknown!", call.=FALSE)
  
  # design characteristics
  desi <- .design.props(d.no)
  #degrees of freedom as expression
  if (robust) {
    dfe  <- parse(text=desi$df2[1],srcfile=NULL)
  } else {
    dfe  <- parse(text=desi$df[1],srcfile=NULL) 
  }

  tval <- qt(1-alpha,eval(dfe))
  s1   <- (log(point)-log(lower))/sqrt(desi$bk/n)/tval
  s2   <- (log(upper)-log(point))/sqrt(desi$bk/n)/tval
  sw   <- 0.5*(s1+s2)
  # both estimates very different?
  if (abs(s1-s2)/sw > 0.1) warning("sw1, sw2 very different. Check input.")
  return(se2CV(sw))
}

# alias to CVfromCI
CI2CV <- function(point, lower, upper, n, design="2x2", alpha=0.05, 
                  robust=FALSE)
{
  if (missing(lower) | missing(upper)) {
    stop("Lower and upper CL must be given!", call.=FALSE) 
  }
  if (missing(point)) point <- sqrt(lower*upper) 
  CVfromCI(point=point, lower, upper, n, design=design, alpha=alpha, 
           robust=robust)
}

