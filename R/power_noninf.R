###############################################################################
# Power & sample size calculations for non-inferiority t-test
# 
# Author: dlabes
###############################################################################

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
    lmarg  <- log(margin)
    diffm  <- log(theta0)
    se     <- CV2se(CV)
  } else {
    if (missing(margin)) margin <- -0.2
    if (missing(theta0)) theta0 <- -0.05
    lmarg  <- margin
    diffm  <- theta0
    se     <- CV
  }
  return(.power.noninf(alpha, lmargin=lmarg, diffm, se, n, df, bk))
}
# --------------------------------------------------------------------------
# internal functions:
# working horse of power function
.power.noninf <- function(alpha=0.025, lmargin, diffm, se, n, df, bk=2)
{
  tval <- qt(1-alpha,df)
  tau  <- abs( (diffm-lmargin)*sqrt(n)/sqrt(bk*se^2) )
  return(1 - pt(tval, df, tau))
}

# start value for sample size search
.sampleN0.noninf <- function(alpha=0.025, targetpower=0.8, lmarg, d0, se, 
                             steps=2, bk=2)
{ 
  n0 <- bk*se^2*(qnorm(targetpower)+ qnorm(1-alpha))^2 / (d0 - lmarg)^2
  n0 <- steps*trunc(n0/steps)
  if (n0<4) n0 <- 4   # minimum sample size
  return(n0)
}
# --------------------------------------------------------------------------
# Sample size estimation for non-inferiority t-test
sampleN.noninf <- function(alpha=0.025, targetpower=0.8, logscale=TRUE, 
                           margin, theta0, CV, design="2x2", robust=FALSE,
                           details=FALSE, print=TRUE, imax=100)
{ 
  #number of the design and check
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design, " unknown!", call.=FALSE)
  
  # design characteristics
  ades   <- .design.props(d.no)
  d.name <- ades$name  # nice name of design
  # get the df for the design as an unevaluated expression (now with n as var)
  if (robust){
    dfe  <- parse(text=ades$df2,srcfile=NULL)
  } else {
    dfe  <- parse(text=ades$df,srcfile=NULL)
  }
  steps  <- ades$steps	# stepsize for sample size search
  bk     <- ades$bk    # get design constant
  
  if (missing(CV)) stop("CV must be given!", call.=FALSE)
  
  if (print) {
    cat("\n+++++++++++ Non-inferiority t-test ++++++++++++\n")
    cat("            Sample size estimation\n")
    cat("-----------------------------------------------\n")
    cat("Study design: ",d.name,"\n")
    if (details) { 
      cat("Design characteristics:\n")
      if (robust & (ades$df2 != ades$df)) {
        cat("df = ",ades$df2," (robust)", sep="") 
      } else cat("df = ",ades$df, sep="")
      cat(", design const. = ", bk, ", step = ", steps,"\n\n",sep="")
    }
  }
  # handle the log transformation
  if (logscale) {
    if (missing(margin)) margin <- 0.8
    if (missing(theta0)) theta0 <- 0.95
    if ( (theta0<=margin) & (margin<1) ) {
      stop("Null ratio ",theta0," must be above margin ",margin,"!", 
          call.=FALSE)
    }
    if ( (theta0>=margin) & (margin>1) ) {
      stop("Null ratio ",theta0," must be below margin ",margin,"!", 
          call.=FALSE)
    }
    lmarg  <- log(margin)
    diffm  <- log(theta0)
    se     <- CV2se(CV)
    if (print) cat("log-transformed data (multiplicative model)\n\n")
  } else {
    if (missing(margin)) margin <- -0.2
    if (missing(theta0)) theta0 <- -0.05
    if ( (theta0<=margin) & (margin<0) ) {
      stop("Null diff. ",theta0," must be above margin ",margin,"!", call.=FALSE)
    }
    if ( (theta0>=margin) & (margin>0) ) {
      stop("Null diff. ",theta0," must be below margin ",margin,"!", call.=FALSE)
    }
    lmarg  <- margin
    diffm  <- theta0
    se     <- CV
    if (print) cat("untransformed data (additive model)\n\n")
  }
  if (print) {
    cat("alpha = ",alpha,", target power = ", targetpower,"\n", sep="")
    cat("Non-inf. margin   = ", margin, "\n", sep="")
    if (logscale) cat("Null (true) ratio = ",theta0,",  CV = ",CV,"\n", sep="")
    else          cat("Null (true) diff. = ",theta0,",  CV = ",CV,"\n", sep="")
  }
  # start value of 'brute force'
  n   <- .sampleN0.noninf(alpha, targetpower, lmarg, d0=diffm, se, steps, bk)
  df  <- eval(dfe)
  pow <- .power.noninf(alpha, lmargin=lmarg, diffm, se, n, df, bk)
  if (details){
    cat("\nSample size search (ntotal)\n")
    cat(" n     power\n")
    # do not print first too high
    # this is for cases with only one step-down and than step up
    if (pow<=targetpower) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  iter <- 0
  while (pow>targetpower){
    if (n<=4) { # min number
      if (details & iter==0) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
      break
    }
    n <- n - steps
    df  <- eval(dfe)
    pow <- .power.noninf(alpha, lmargin=lmarg, diffm, se, n, df, bk)
    iter <- iter+1
    if (details) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    if (iter>imax) break  
  }
  while(pow<targetpower){
    n <- n+steps
    df  <- eval(dfe)
    pow <- .power.noninf(alpha, lmargin=lmarg, diffm, se, n, df, bk)
    if (details) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    iter <- iter +1
    if (iter>imax) break  
  }
  if (pow<targetpower) {
    n <- NA
    if (details) cat("Sample size search failed!\n")
  }
  if (print && !details) {
    cat("\nSample size (total)\n")
    cat(" n     power\n")
    cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  # return value: a data.frame
  res <- data.frame(design=design, alpha=alpha, CV=CV, theta0=theta0, 
         margin=margin, n=n, power=pow, targetpower=targetpower)
  names(res) <- c("Design","alpha","CV","theta0","Margin", "Sample size", 
                  "Achieved power", "Target power")
  
  if (print) return(invisible(res)) else return(res)
}
