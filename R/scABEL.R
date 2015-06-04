# helper function for the scABEL regulatory settings 
reg_const <- function(regulator)
{
  if (regulator=="FDA"){
    r <- list(CVswitch=0.3, r_const=log(1.25)/0.25, CVcap=Inf)
  } 
  if (regulator=="ANVISA"){
    # same regulatory const. as EMA but
    # switch to widened limits if CVRef>=40%
    r <- list(CVswitch=0.4, r_const=log(1.25)/CV2se(0.3), CVcap=0.5)
  }
  if (regulator=="EMA"){
    # r_const taken literally from BE guideline
    r <- list(CVswitch=0.3, r_const=0.76, CVcap=0.5)
  }
  r
}  
# -------------------------------------------------------------------------
# function to calculate the widened ABE limits
scABEL <- function(CV, regulator=c("EMA", "ANVISA", "FDA", "USER"), 
                   r_const, CVswitch, CVcap)
{
  regulator <- match.arg(regulator)
  if (regulator=="USER"){
    if (missing(CVswitch) | missing(CVcap) | missing(r_const)){
      stop("r_const, CVswitch, CVcap must be given.")
    }
  } else {
    rc    <- reg_const(regulator)
    CVcap <- rc$CVcap
    CVswitch <- rc$CVswitch
    r_const  <- rc$r_const
  }

  ret <- ifelse(CV<=CVswitch, 1.25, exp(r_const*CV2se(CV)))
  ret <- ifelse(CV>CVcap, exp(r_const*CV2se(CVcap)), ret)
  if (length(CV)>1){
    ret <- cbind(1/ret, ret)
    colnames(ret) <- c("lower", "upper")
  } else {
    ret <- c(1/ret, ret)
    names(ret) <- c("lower", "upper")
  }
  ret
}