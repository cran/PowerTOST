#---------------------------------------------------------------------------
# Sample size for partial and full replicate design and scaled ABE 
# via simulated (empirical) power
# 
# Author: dlabes
#---------------------------------------------------------------------------

sampleN.scABEL <- function(alpha=0.05, targetpower=0.8, theta0, theta1, 
                           theta2, CV, design=c("2x3x3", "2x2x4"), 
                           regulator=c("EMA", "FDA"), nsims=1E5,
                           nstart, print=TRUE, details=TRUE, setseed=TRUE)
{
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta0)) theta0 <- 0.95
  if (missing(theta2)) theta2=1/theta1
  if ( (theta0<=theta1) | (theta0>=theta2) ) {
    stop("Null ratio ",theta0," not between margins ",theta1," / ",theta2,"!", 
         call.=FALSE)
  }
  if (missing(CV)) stop("CV must be given!", call.=FALSE)
  
  # for later enhancement taking into account the 
  # subject-by-formulation interaction
  # can we incorporate this? EMA method doesn't have such a term.
  s2D  <- 0  
  CVwT <- CV[1]
  # should we allow different variabilities in the EMA method?
  if (length(CV)==2) CVwR <- CV[2] else CVwR <- CVwT
  s2wT <- log(1.0 + CVwT^2)
  s2wR <- log(1.0 + CVwR^2)
  
  regulator <- match.arg(regulator)
  if (regulator=="FDA"){
    CVcap    <- Inf
    CVswitch <- 0.3
    r_const  <- log(1.25)/0.25
  } else {
    # regulatory settings for EMA
    CVcap    <- 0.5
    CVswitch <- 0.3
    r_const  <- 0.760
  }
  
  # check design
  design <- match.arg(design)
  # we are treating only balanced designs
  # thus we use here bk - design constant for ntotal
  # expressions for the df's
  if (design=="2x3x3") {
    bk <- 1.5; seqs <- 3
    dfe   <- parse(text="2*n-3", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    #sd2  <- s2D + (s2wT + s2wR)/2 # used in v1.1-00 - v1.1-02
    # simulations with s2D=0 show:
    Emse  <- (s2wT + 2.0*s2wR)/3
    cvec  <- c(1, 2) # for sim of mses from s2wT and s2wR
  }
  if (design=="2x2x4") {
    bk <- 1.0; seqs <- 2
    # only EMA settings
    dfe   <- parse(text="3*n-4", srcfile=NULL)
    dfRRe <- parse(text="n-2", srcfile=NULL)
    # sd^2 (variance) of the differences T-R from their components
    Emse  <- (s2wT + s2wR)/2
    cvec  <- c(1, 1)
  }
  mlog <- log(theta0)
  
  if (print){
    cat("\n+++++++++++ scaled (widened) ABEL +++++++++++\n")
    cat("            Sample size estimation\n")
    cat("---------------------------------------------\n")
    cat("Study design: ",design,"\n")
    cat("log-transformed data (multiplicative model)\n")
    cat(nsims,"studies simulated.\n\n")
    cat("alpha  = ",alpha,", target power = ", targetpower,"\n", sep="")
    cat("CVw(T) = ",CVwT,"; CVw(R) = ",CVwR,"\n", sep="")
    cat("Null (true) ratio = ",theta0,"\n", sep="")
    cat("ABE limits / PE constraints =",theta1,"...", theta2,"\n")
    cat("Regulatory settings:",regulator,"\n")
    if (details) { 
      cat("- CVswitch = ", CVswitch)
      if (is.finite(CVcap)){
        cat(", cap on ABEL if CV > ", CVcap,"\n",sep="")
      } else {
        cat(", no cap on ABEL\n",sep="")
      }  
      cat("- Regulatory constant =",r_const,"\n")
    }     
  }
  
  # -----------------------------------------------------------------
  # nstart? from sampleN0 with widened limits
  # does'nt fit really good if theta0>=1.2! ways out?
  ltheta1 <- -sqrt(s2wR)*r_const
  ltheta2 <- -ltheta1
  if (CVwR <= CVswitch){
    ltheta1 <- log(theta1)
    ltheta2 <- log(theta2)
  }
  if (CVwR > CVcap){
    ltheta1 <- -sqrt(log(1.0 + CVcap^2))*r_const
    ltheta2 <- -ltheta1
  }
  if (missing(nstart)){
    n <- .sampleN0(alpha=alpha, targetpower, ltheta1, ltheta2, diffm=mlog, 
                   se=sqrt(Emse), steps=seqs, bk=bk)
  } else n <- seqs*round(nstart/seqs)           
  # iterate until pwr>=targetpower
  # we are simulating for balanced designs
  C2 <- bk/n
  # sd of the sample mean T-R (point estimator)
  sdm  <- sqrt(Emse*C2)
  df   <- eval(dfe)
  dfRR <- eval(dfRRe)
  
  if(setseed) set.seed(123456)
  p <- .power.scABEL(mlog, sdm, C2, Emse, cvec, df, s2wR, dfRR, s2wT,
                     nsims, CVswitch, r_const, CVcap, 
                     ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha=alpha)
  pwr <- as.numeric(p["BE"]);
  
  if (details) {
    cat("\nSample size search\n")
    cat(" n     power\n")
    # do not print first too high
    # this is for cases with only one step-down and than step up
    if (pwr<=targetpower) cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
  }
  iter <- 0; imax=100
  nmin <- 6
  # iter>100 is emergency brake
  # --- loop until power <= target power, step-down
  while (pwr>targetpower) {
    if (n<=nmin) { 
      if (details & iter==0) cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
      break
    }
    n  <- n-seqs     # step down if start power is to high
    iter <- iter + 1
# rising the stepsize does not function stepsize too big in some cases    
#    if (abs(pwr-targetpower)>0.03) n  <- n-seqs
    C2 <- bk/n
    # sd of the sample mean T-R (point estimator)
    sdm  <- sqrt(Emse*C2)
    df   <- eval(dfe)
    dfRR <- eval(dfRRe)
    
    if(setseed) set.seed(123456)
    p <- .power.scABEL(mlog, sdm, C2, Emse, cvec, df, s2wR, dfRR, s2wT,
                       nsims, CVswitch, r_const, CVcap, 
                       ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha=alpha)
    pwr <- as.numeric(p["BE"]);
    
    if (details) cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
    if (iter>imax) break  
    # loop results in n with power too low
    # must step one up again. is done in the next loop
  }
  while (pwr<targetpower) {
    n    <- n+seqs   # step-up
#    does not function stepsize too big in some cases    
#    if (abs(pwr-targetpower)>0.03) n  <- n+seqs
    iter <- iter+1
    C2 <- bk/n
    sdm  <- sqrt(Emse*C2)
    df   <- eval(dfe)
    dfRR <- eval(dfRRe)
    
    if(setseed) set.seed(123456)
    p <- .power.scABEL(mlog, sdm, C2, Emse, cvec, df, s2wR, dfRR, s2wT, 
                       nsims, CVswitch, r_const, CVcap, 
                       ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha=alpha)
    pwr <- as.numeric(p["BE"]);
    
    if (details) cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
    if (iter>imax) break 
  }
  
  if (pwr<targetpower) {
    n <- NA
    if (details) cat("Sample size search failed!\n")
  }
  
  if (print && !details) {
    cat("\nSample size\n")
    cat(" n     power\n")
    cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
    if (is.na(n)) cat("Sample size search failed!\n")
  }
  if (print) cat("\n")
  
  if (print) return(invisible(n)) 
  else return(n)
  
} # end function


