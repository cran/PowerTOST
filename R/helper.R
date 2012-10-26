#------------------------------------------------------------------------------
# Author: dlabes
#------------------------------------------------------------------------------

# helper function to calculate std err from CV of lognormal data
CV2se <- function(CV) sqrt(log(1.0 + CV^2))
# reverse to CV2se helper function
se2CV <- function(se) sqrt(exp(se*se)-1)
# helper function to calculate mse from CV of lognormal data
CV2mse <- function(CV) log(1.0 + CV^2)
mse2CV <- function(mse) sqrt(exp(mse)-1)

# function to calculate confidence limits of given CV 
CVCL <- function(CV, df, side=c("upper", "lower","2-sided"), alpha=0.05)
{
  ssintra <- log(1.0 + CV^2)*df
  side    <- match.arg(side)
  if(side=="upper") return(sqrt(exp(ssintra/qchisq(alpha,df))-1))
  if(side=="lower") return(sqrt(exp(ssintra/qchisq(1-alpha,df))-1))
  if(side=="2-sided") {
    limits <- c(ssintra/qchisq(1-alpha/2,df),ssintra/qchisq(alpha/2,df))
    return(sqrt(exp(limits)-1))
  }  
}
