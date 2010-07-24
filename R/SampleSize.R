# Author: dlabes
# ----- helper functions for sampleN.TOST -------------------------------------
# Sample size for a desired power, large sample approx.
# bk = design constant, see known.designs()
.sampleN0 <- function(alpha=0.05, targetpower=0.8, ltheta1, ltheta2, diffm, se, 
		                  steps=2, bk=2)
{
  z1 <- qnorm(1-alpha)
  if (abs(diffm)>0.0001) z2 <- qnorm(targetpower) 
    else z2 <- qnorm(1-(1-targetpower)/2)
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
# leave upper BE margin (theta2) empty and the function will use -lower
# in case of additive model or 1/lower if logscale=TRUE
sampleN.TOST <- function(alpha=0.05, targetpower=0.8, logscale=TRUE, 
                         theta1=0.8, theta2, diff=0.95, CV, 
                         design="2x2", exact=TRUE, print=TRUE, details=FALSE)
{
  #number of the design and check
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Err: design not known", call.=FALSE)
  
  # design characteristics
  ades   <- .design.props(d.no)
  d.name <- ades$name  # nice name of design
  # get the df for the design as an unevaluated expression (now with n as var)
  dfe    <- parse(text=ades$df,srcfile=NULL) 
  steps  <- ades$steps	# stepsize for sample size search
  bk     <- ades$bk # get design constant
  
  if (missing(CV)) stop("Err: CV must be given!", call.=FALSE)
  
  # print the configuration:
  if (print) {
    cat("\n+++++++++ Equivalence test - TOST +++++++++\n")
    cat("          Sample size estimation\n")
    cat("-------------------------------------------\n")
    cat("Study design: ",d.name,"\n")
    if (details) { 
      cat("Design characteristics:\n")
      cat("df = ",ades$df,", design const. = ",bk,
          ", step = ",steps,"\n\n",sep="")
    }     
  }
  # handle the log transformation
  if (logscale) {
    if (missing(theta2)) theta2=1/theta1
    if ( (diff<=theta1) | (diff>=theta2) ) {
      stop("Err: ratio not between margins!", call.=FALSE)
    }
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
    diffm   <- log(diff)
    se      <- CV2se(CV)
    if (print) cat("log-transformed data (multiplicative model)\n\n")
  } else {
    if (missing(theta2)) theta2=-theta1
    if ( (diff<=theta1) | (diff>=theta2) ) {
      stop("Err: diff not between margins!", call.=FALSE)
    }
    ltheta1 <- theta1
    ltheta2 <- theta2
    diffm   <- diff
    se      <- CV
    if (print) cat("untransformed data (additive model)\n\n")
  }
  
  if (print) {
    cat("alpha = ",alpha,", target power = ", targetpower,"\n", sep="")
    cat("BE margins        =",theta1,"...", theta2,"\n")
    if (logscale) cat("Null (true) ratio = ",diff,",  CV = ",CV,"\n", sep="")
    else          cat("Null (true) diff. = ",diff,",  CV = ",CV,"\n", sep="")
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
    cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  # --- loop until power >= target power
  iter <- 0
  # iter>50 is emergency brake
  # this is eventually not necessary, depends on quality of sampleN0
  # in experimentation I have seen max of six steps
  if (pow<targetpower) {
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
  } else {
    while (pow>targetpower) {
      if (n<=4) break     # min. sample size
      n    <- n-steps     # step down if start power is to high
      iter <- iter+1
      df   <- eval(dfe)
      if (exact)
        pow  <- .power.TOST(alpha, ltheta1,ltheta2,diffm, se, n, df, bk) else
        pow  <- .approx.power.TOST(alpha,ltheta1,ltheta2, diffm, se, n, df, bk)
      if (details) cat( n," ", formatC(pow, digits=6),"\n")
      if (iter>50) break  
    }
    if ( pow<targetpower ) {
      n   <- n+steps     #step up once to have n with pow>=target
      df  <- eval(dfe)
      if (exact)
        pow  <- .power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk) else
        pow  <- .approx.power.TOST(alpha, ltheta1, ltheta2, diffm, se, n, df, bk)
    }
  }
  if (print && !details) {
    cat("\nSample size\n")
    if (d.no == 0) cat("(n is sample size per group)\n") #parallel group design
    cat(" n     power\n")
    cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  if (details && print) {
    if (exact) cat("\nExact power calculation with\nOwen's Q functions.\n")
    else 
      cat("\nApproximate power calculation with\nnon-central t-distribution.\n")
  }
  # always print if approx.
  if (print && !exact) 
    cat("\nApproximate power calculation with\nnon-central t-distribution.\n")
  if (print) cat("\n")
  
  #return results as data.frame
  res <- data.frame(design=design, alpha=alpha, CV=CV, theta0=diff, 
      theta1=theta1, theta2=theta2, n=n, power=pow, 
      targetpower=targetpower)
  names(res) <-c("Design","alpha","CV","theta0","theta1","theta2",
      "Sample size", "Achieved power", "Target power")
  
  if (print) return(invisible(res)) 
  else return(res)
  
}
#------------------------------------------------------------------------------
# sample size start for Joulious "expected" power
.expsampleN0 <- function(alpha=0.05, targetpower=0.8, ltheta1, ltheta2, diffm, 
                         se, dfse, steps=2, bk=2)
{
  Z1 <- qnorm(1-alpha)
  if (abs(diffm)>0.0001) tinv <- qt(targetpower, dfse, Z1)  else
    tinv <- qt(1-(1-targetpower)/2, dfse, Z1) 
  
  # is factor 2 in julious = bk?
  n01  <- bk*(se*tinv/(ltheta1-diffm))^2
  n02  <- bk*(se*tinv/(ltheta2-diffm))^2
  # print(n01);print(n02)
  n0 <- ceiling(max(n01,n02))
  #make an even multiple of step (=2 in case of 2x2 cross-over)
  n0 <- steps*trunc(n0/steps)
  if (n0<4) n0 <- 4   # minimum sample size
  
  return(n0)
}
#------------------------------------------------------------------------------
# Sample size for a desired "expected" power according to Julious: 
# see known.designs() for covered experimental designs
# Only for log-transformed data
# leave upper BE margin (ltheta2) empty and the function will use 1/lower
# CV and dfCV can be vectors, if then a pooled CV, df will be calculated
expsampleN.TOST <- function(alpha=0.05, targetpower=0.8, 
                            theta1=0.8, theta2, diff=0.95, CV, dfCV, 
                            alpha2=0.05,
                            design="2x2", print=TRUE, details=FALSE)
{
  #number of the design and check if design is implemented
  d.no <- .design.no(design)
  if (is.na(d.no)) stop("Err: design not known", call.=FALSE)
  
  # design characteristics
  ades   <-.design.props(d.no)
  d.name <- ades$name  # nice name of design
  # get the df for the design as an unevaluated expression
  dfe    <- parse(text=ades$df,srcfile=NULL) 
  steps  <- ades$steps	# stepsize for sample size search
  bk     <- ades$bk # get design constant
  
  if (missing(CV) | missing(dfCV)) {
    stop("Err: CV  and df must be given!", call.=FALSE)
  }
  
  # calculate pooled data if CV and dfCV are vectors
  if (length(CV)>1){
    if (length(dfCV)!=length(CV)) {
      stop("Err: CV and df must have equal number of entries!", call.=FALSE)
    }
    dfse <- sum(dfCV)
    CVp  <- CV2se(CV)^2 #need s-squared
    CVp  <- CVp * dfCV
    CVp  <- sum(CVp)/dfse
    CVp  <- sqrt(CVp)
    CVp  <- se2CV(CVp)
  } else {
    dfse <- dfCV
    CVp  <- CV
  }
  # print the configuration:
  if (print) {
    cat("\n++++++++ Equivalence test - TOST ++++++++\n")
    cat("   Sample size est. with uncertain CV\n")
    cat("-----------------------------------------\n")
    cat("Study design: ",d.name,"\n")
    if (details) { 
      cat("Design characteristics:\n")
      cat("df = ",ades$df,", design const. = ",bk,", step = ",steps,
          "\n\n",sep="")
    }     
  }
  if (missing(theta2)) theta2=1/theta1
  if ( (diff<=theta1) | (diff>=theta2) ) {
      stop("Err: ratio not between margins!", call.=FALSE)
    }
  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  diffm   <- log(diff)
  se      <- CV2se(CVp)
  
  if (print) {
    cat("log-transformed data (multiplicative model)\n\n")
    cat("alpha = ",alpha,", target power = ", targetpower,"\n", sep="")
    cat("BE margins        =",theta1,"...", theta2,"\n")
    cat("Null (true) ratio = ",diff,"\n", sep="")
    # can use lower.tail=FALSE and 1-0.05 in qchisq, H. Schütz in his lectures
    seupper <- dfse*se^2/qchisq(alpha2, dfse)
    seupper <- sqrt(seupper)
    if (length(CV)>1){
      cat("Variability data\n")
      print(data.frame(CV=CV,df=dfCV), row.names = FALSE)
      cat("CV(pooled)         = ", CVp, " with ", dfse," df\n", sep="")
    } else {
      cat("CV                 = ", CVp, " with ", dfse," df\n", sep="")
    }   
    cat("one-sided upper CL = ",se2CV(seupper)," (level = ",
        100*(1-alpha2),"%)\n",sep="")
  }
  
  #start value from large sample approx. 
  n   <- .expsampleN0(alpha, targetpower, ltheta1, ltheta2, diffm, 
                      se, dfse, steps, bk)
  df  <- eval(dfe)
  
  pow <- .exppower.TOST(alpha, ltheta1, ltheta2, diffm, se, dfse, n, df, bk) 
  if (details) {
    cat("\nSample size search\n")
    if (d.no == 0) cat("(n is sample size per group)\n") #parallel group design
    cat(" n    exp. power\n")
    cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  # --- loop until power >= target power
  iter <- 0
  # iter>50 is emergency brake
  # this is eventually not necessary, depends on quality of sampleN0
  # in experimentation I have seen max of six steps
  if (pow<targetpower) {
    while (pow<targetpower) {
      n    <- n+steps
      iter <- iter+1
      df   <- eval(dfe)
      pow <- .exppower.TOST(alpha, ltheta1, ltheta2, diffm, se, dfse,
                             n, df, bk) 
      if (details) cat( n," ", formatC(pow, digits=6, format="f"),"\n")
      if (iter>50) break 
    }
  } else {
    while (pow>targetpower) {
      if (n<=4) break     # min. sample size
      n    <- n-steps     # step down if start power is to high
      iter <- iter+1
      df   <- eval(dfe)
      pow  <- .exppower.TOST(alpha, ltheta1, ltheta2, diffm, se, dfse, 
                              n, df, bk) 
      if (details) cat( n," ", formatC(pow, digits=6),"\n")
      if (iter>50) break  
    }
    if ( pow<targetpower ) {
      n    <- n+steps     #step up once to have n with pow>=target
      df   <- eval(dfe)
      pow  <- .exppower.TOST(alpha, ltheta1, ltheta2, diffm, se, dfse, 
                              n, df, bk) 
    }
  }
  if (print && !details) {
    cat("\nSample size\n")
    if (d.no == 0) cat("(n is sample size per group)\n") #parallel group design
    cat(" n    exp. power\n")
    cat( n," ", formatC(pow, digits=6, format="f"),"\n")
  }
  
  if (print) cat("\n")
  
  #return results as data.frame
  res <- data.frame(design=design, alpha=alpha, CV=CV, dfCV=dfse, theta0=diff, 
                    theta1=theta1, theta2=theta2, n=n, power=pow, 
                    targetpower=targetpower)
  names(res) <-c("Design","alpha","CV","df of CV","theta0","theta1","theta2",
                 "Sample size", "Achieved power", "Target power")
  
  if (print) return(invisible(res)) 
  else return(res)
  
}

