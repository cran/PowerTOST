# Author: dlabes
# ----- helper functions for sampleN.TOST -------------------------------------
# Sample size for a desired power, large sample approx.
# bk = design constant, see known.designs()
.sampleN0 <- function(alpha=0.05, targetpower=0.8, ltheta1, ltheta2, diffm, se, 
		                  steps=2, bk=2)
{
  z1 <- qnorm(1-alpha)
  if (abs(diffm)>0.0001) z2 <- qnorm(targetpower) else
     z2 <- qnorm(1-(1-targetpower)/2)
  n01<-(bk/2)*( (z1+z2)*(se*sqrt(2)/(diffm-ltheta1)) )^2;
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
sampleN.TOST <- function(alpha=0.05, targetpower=0.8, logscale=TRUE, 
                         theta1, theta2, theta0, CV, design="2x2",
                         exact=TRUE, print=TRUE, details=FALSE)
{
  #number of the design and check
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Design not known!", call.=FALSE)
  
  # design characteristics
  ades   <- .design.props(d.no)
  d.name <- ades$name  # nice name of design
  # get the df for the design as an unevaluated expression (now with n as var)
  dfe    <- parse(text=ades$df,srcfile=NULL) 
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
      cat("df = ",ades$df,", design const. = ",bk,
          ", step = ",steps,"\n\n",sep="")
    }     
  }
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
  if (exact) 
    pow <- .power.TOST(alpha, ltheta1, ltheta2, diffm, se, n=n, df, bk) else
    pow <- .approx.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n=n, df, bk)
  if (details) {
    cat("\nSample size search\n")
    if (d.no == 0) cat("(n is sample size per group)\n") #parallel group design
    cat(" n     power\n")
    # do not print first too high
    if (pow<=targetpower) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  iter <- 0
  # iter>50 is emergency brake
  # this is eventually not necessary, depends on quality of sampleN0
  # in experimentation I have seen max of six steps
  while (pow>targetpower) {
    if (n<=4) { # min number
      if (print & iter==0) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
      break
    }
    n    <- n-steps     # step down if start power is to high
    iter <- iter+1
    df   <- eval(dfe)
    if (exact)
      pow  <- .power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk) else
      pow  <- .approx.power.TOST(alpha,ltheta1,ltheta2,diffm,se,n,df,bk)
    # do not print first step down
    if (details) cat( n," ", formatC(pow, digits=6),"\n")
    if (iter>50) break  
    # loop results in n with power too low
    # must step one up again
  }
  # --- loop until power >= target power
  while (pow<targetpower) {
    n    <- n+steps
    iter <- iter+1
    df   <- eval(dfe)
    if (exact)
      pow  <- .power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk) else
      pow  <- .approx.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk)
    if (details) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
    if (iter>50) break 
  }
  if (print && !details) {
    cat("\nSample size\n")
    if (d.no == 0) cat("(n is sample size per group)\n") #parallel group design
    cat(" n     power\n")
    cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  if (details && print) {
    if (exact) cat("\nExact power calculation with\nOwen's Q functions.\n")
  }
  # always print if approx.
  if (print && !exact) 
    cat("\nApproximate power calculation with\nnon-central t-distribution.\n")
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
