# Author: dlabes
# ----- helper functions for sampleN.TOST -------------------------------------
# Sample size for a desired power, large sample approx.
# bk = design constant, see known.designs()
.sampleN0 <- function(alpha=0.05, targetpower=0.8, ltheta1, ltheta2, diffm, se, 
		                  steps=2, bk=2)
{
  z1 <- qnorm(1-alpha)
  # value 0.04 corresponds roughly to log(0.96)
  # with lower values there are many steps around between 0.95 and 1
  if (abs(diffm)>0.04) z2 <- qnorm(targetpower) else {
    z2 <- qnorm(1-(1-targetpower)/2) # diffm ~0 (log: theta0=1)
    diffm <- 0
  }
  n01<-(bk/2)*((z1+z2)*(se*sqrt(2)/(diffm-ltheta1)))^2;
  n02<-(bk/2)*((z1+z2)*(se*sqrt(2)/(diffm-ltheta2)))^2;
  n0 <- ceiling(max(n01,n02))
  #make an even multiple of step (=2 in case of 2x2 cross-over)
  n0 <- steps*trunc(n0/steps)
  if (n0<4) n0 <- 4   # minimum sample size
	
  return(n0)
}

#------------------------------------------------------------------------------
# Sample size for a desired power: 
# see known.designs() for covered experimental designs
# theta1 if empty is set to 0.8 or -0.2 depending on logscale
# diff if empty is set to 0.95 or 0.05 depending on logscale
# leave upper BE margin (theta2) empty and the function will use -lower
# in case of additive model or 1/lower if logscale=TRUE
sampleN.TOST <- function(alpha=0.05, targetpower=0.8, logscale=TRUE, theta0, 
                         theta1, theta2, CV, design="2x2", method="exact",
                         robust=FALSE, print=TRUE, details=FALSE, imax=100)
{
  #number of the design and check
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design ",design, " unknown!", call.=FALSE)
  
  # design characteristics
  ades   <- .design.props(d.no)
  d.name <- ades$name  # nice name of design
  # get the df for the design as an unevaluated expression (now with n as var)
  if (robust){
    dfe    <- parse(text=ades$df2,srcfile=NULL)
  } else {
    dfe    <- parse(text=ades$df,srcfile=NULL)
  }
  steps  <- ades$steps	# stepsize for sample size search
  bk     <- ades$bk # get design constant
  
  if (missing(CV)) stop("CV must be given!", call.=FALSE)
  
  # print the configuration:
  if (print) {
    cat("\n+++++++++++ Equivalence test - TOST +++++++++++\n")
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
  
  # regularize the method giving
  method <- .powerMethod(method)
  
  # handle the log transformation
  if (logscale) {
    if (missing(theta1) & missing(theta2)) theta1 <- 0.8
    if (missing(theta0)) theta0 <- 0.95
    if (missing(theta2)) theta2=1/theta1
    if ( (theta0<=theta1) | (theta0>=theta2) ) {
      stop("Null ratio ",theta0," not between margins ",theta1," / ",theta2,"!", 
           call.=FALSE)
     }
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    diffm   <- log(theta0)
    se      <- CV2se(CV)
    if (print) cat("log-transformed data (multiplicative model)\n\n")
  } else {
    if (missing(theta1) & missing(theta2)) theta1 <- -0.2
    if (missing(theta0)) theta0 <- 0.05
    if (missing(theta2)) theta2=-theta1
    if ( (theta0<=theta1) | (theta0>=theta2) ) {
      stop("Null diff. ",theta0," not between margins ",theta1," / ",theta2,"!", 
           call.=FALSE)
    }
    ltheta1 <- theta1
    ltheta2 <- theta2
    diffm   <- theta0
    se      <- CV
    if (print) cat("untransformed data (additive model)\n\n")
  }
  
  if (print) {
    cat("alpha = ",alpha,", target power = ", targetpower,"\n", sep="")
    cat("BE margins        =",theta1,"...", theta2,"\n")
    if (logscale) cat("Null (true) ratio = ",theta0,",  CV = ",CV,"\n", sep="")
    else          cat("Null (true) diff. = ",theta0,",  CV = ",CV,"\n", sep="")
  }
  
  # start value from large sample approx. (hidden func.)
  n  <- .sampleN0(alpha, targetpower, ltheta1, ltheta2, diffm, se, steps, bk)
  df <- eval(dfe)
  pow <- .calc.power(alpha, ltheta1, ltheta2, diffm, se, n, df, bk, method)

  if (details) {
    cat("\nSample size search (ntotal)\n")
    # parallel group design is now in terms of ntotal
    #if (d.no == 0) cat("(n is sample size per group)\n") #parallel group design
    cat(" n     power\n")
    # do not print first too high
    # this is for cases with only one step-down and than step up
    if (pow<=targetpower) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  iter <- 0
  # iter>100 is emergency brake
  # this is eventually not necessary, depends on quality of sampleN0
  # in experimentation I have seen max of six steps
  # reformulation with only one loop does not shorten the code considerable
  # --- loop until power <= target power, step-down
  while (pow>targetpower) {
    if (n<=4) { # min number
      if (details & iter==0) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
      break
    }
    n    <- n-steps     # step down if start power is to high
    iter <- iter+1
    df   <- eval(dfe)
    pow <- .calc.power(alpha, ltheta1, ltheta2, diffm, se, n, df, bk, method)
    
    # do not print first step down
    if (details) cat( n," ", formatC(pow, digits=6),"\n")
    if (iter>imax) break  
    # loop results in n with power too low
    # must step one up again. is done in the next loop
  }
  # --- loop until power >= target power
  while (pow<targetpower) {
    n    <- n+steps
    iter <- iter+1
    df   <- eval(dfe)
    pow <- .calc.power(alpha, ltheta1, ltheta2, diffm, se, n, df, bk, method)
    if (details) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    if (iter>imax) break 
  }

  if (pow<targetpower) {
    n <- NA
    if (details) cat("Sample size search failed!\n")
  }
  
  if (print && !details) {
    cat("\nSample size (total)\n")
    #if (d.no == 0) cat("(n is sample size per group)\n") #parallel group design
    cat(" n     power\n")
    cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  if (details && print) {
    if (method=="exact") 
      cat("\nExact power calculation with\nOwen's Q functions.\n")
  }
  # always print if approx.
  if (print && method!="exact"){
    approx <- switch(
      method,
      noncentral="Approximate power calculation with\nnon-central t-distribution.",
      shifted="Approximate power calculation with\nshifted central t-distribution."
      )
    cat("\n",approx,"\n",sep="")
  } 
  if (print) cat("\n")
  
  #return results as data.frame
  res <- data.frame(design=design, alpha=alpha, CV=CV, theta0=theta0, 
                    theta1=theta1, theta2=theta2, n=n, power=pow, 
                    targetpower=targetpower)
  names(res) <-c("Design","alpha","CV","theta0","theta1","theta2",
                 "Sample size", "Achieved power", "Target power")
  
  if (print) return(invisible(res)) 
  else return(res)
  
}
