#---------------------------------------------------------------------------
# Sample size for partial and full replicate design and scaled ABE 
# via simulated (empirical) power
# 
# Author: dlabes
#---------------------------------------------------------------------------

sampleN.scABEL <- function(alpha=0.05, targetpower=0.8, theta0, theta1, 
                           theta2, CV, design=c("2x3x3", "2x2x4"), 
                           regulator=c("EMA", "FDA"), nsims=1E6,
                           nstart, print=TRUE, details=TRUE)
{
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta0)) theta0 <- 0.95
  if (missing(theta2)) theta2=1/theta1
  if ( (theta0<=theta1) | (theta0>=theta2) ) {
    stop("Null ratio ",theta0," not between margins ",theta1," / ",theta2,"!", 
         call.=FALSE)
  }
  if (missing(CV)) stop("CV must be given!", call.=FALSE)
  
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
    # in case of the FDA we are using the 'robust' df's
    # due to the fact that the described analysis in the
    # progesterone guidance is based in the intrasubject contrasts
    # T-R and R-R with df=n-seqs
    if (regulator=="FDA"){
      dfe   <- parse(text="n-3", srcfile=NULL)
      dfRRe <- parse(text="n-3", srcfile=NULL)
    } else {
      dfe   <- parse(text="2*n-3", srcfile=NULL)
      dfRRe <- parse(text="n-2", srcfile=NULL)
    }
  }
  if (design=="2x2x4") {
    bk <- 1.0; seqs <- 2
    if (regulator=="FDA"){
      dfe   <- parse(text="n-2", srcfile=NULL)
      dfRRe <- parse(text="n-2", srcfile=NULL)
    } else {
      # EMA settings
      dfe   <- parse(text="3*n-4", srcfile=NULL)
      dfRRe <- parse(text="n-3", srcfile=NULL)
    }
  }
  # for later enhancement taking into account the 
  # subject-by-formulation interaction
  sD2  <- 0  
  sWT2 <- log(1.0 + CV[1]^2)
  if (length(CV)==2) sWR2 <- log(1.0 + CV[2]^2) else sWR2 <-sWT2
  CVwR <- sqrt(exp(sWR2)-1)
  CVwT <- sqrt(exp(sWT2)-1)
  # sd^2 of the differences T-R from their components
  sd2  <- (sD2 + (sWT2 + sWR2)/2) # is this correct for partial replicate?
  mlog <- log(theta0)
  
  if (print){
    cat("\n+++++++++++ scaled (widened) ABEL +++++++++++\n")
    cat("            Sample size estimation\n")
    cat("---------------------------------------------\n")
    cat("Study design: ",design,"\n")
    cat("log-transformed data (multiplicative model)\n\n")
    cat("alpha = ",alpha,", target power = ", targetpower,"\n", sep="")
    cat("CVw(T) = ",CVwT,"; CVw(R) = ",CVwR,"\n", sep="")
    cat("Null (true) ratio = ",theta0,"\n", sep="")
    cat("PE constraints    =",theta1,"...", theta2,"\n")
    cat("Regulatory body:",regulator,"\n")
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
  ltheta1 <- -sqrt(sWR2)*r_const
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
                   se=sqrt(sd2), steps=seqs, bk=bk)
  } else n <- nstart           
  # next is an empirical observation for the EMA settings
  # both together give n=2*n in case of CVwR>0.5
  # TODO: check this if CVwR != CVwT
#  if (targetpower<0.9){
#    if ((theta0>=1.2 | theta0<=0.833) & CVwR>=0.45) n <- 1.5*n
#    if ((theta0>=1.2 | theta0<=0.833) & CVwR>=0.5)  n <- 1.3333*n
#  } else {
#    if ((theta0>=1.2 | theta0<=0.833) & CVwR>=0.45) n <- 1.4*n
#    if ((theta0>=1.2 | theta0<=0.833) & CVwR>=0.5)  n <- 1.8*n
#  } 
#  n <- seqs*round(n/seqs, 0)
  # iterate until pwr>=targetpower
  # we are simulating for balanced designs
  fact <- bk/n
  # sd of the sample mean T-R (point estimator)
  sdm  <- sqrt(sd2*fact)
  df   <- eval(dfe)
  dfRR <- eval(dfRRe)
  p <- .power.scABEL(mlog, sdm, fact, sd2, df, sWR2, dfRR, nsims, 
                     ln_lBEL=log(theta1),ln_uBEL=log(theta2), 
                     CVswitch, r_const, CVcap, alpha=alpha)
  pwr <- as.numeric(p["BE"]);
  
  if (details) {
    cat("\nSample size search\n")
    # parallel group design is now in terms of ntotal
    #if (d.no == 0) cat("(n is sample size per group)\n") #parallel group design
    cat(" n     power\n")
    # do not print first too high
    # this is for cases with only one step-down and than step up
    if (pwr<=targetpower) cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
  }
  iter <- 0; imax=100
  nmin <- 6
  # iter>100 is emergency brake
  # this is eventually not necessary, depends on quality of sampleN0
  # in experimentation I have seen max of six steps
  # reformulation with only one loop does not shorten the code considerable
  # --- loop until power <= target power, step-down
  while (pwr>targetpower) {
    if (n<=nmin) { 
      if (details & iter==0) cat( n," ", formatC(pwr, digits=6, format="f"),"\n")
      break
    }
    n  <- n-seqs     # step down if start power is to high
#    does not function stepsize too big in some cases    
#    if (abs(pwr-targetpower)>0.03) n  <- n-seqs
    fact <- bk/n
    # sd of the sample mean T-R (point estimator)
    sdm  <- sqrt(sd2*fact)
    df   <- eval(dfe)
    dfRR <- eval(dfRRe)
    
    p <- .power.scABEL(mlog, sdm, fact, sd2, df, sWR2, dfRR, nsims, 
                       ln_lBEL=log(theta1),ln_uBEL=log(theta2), 
                       CVswitch, r_const, CVcap, alpha=alpha)
    pwr <- as.numeric(p["BE"]);
    
    # do not print first step down
    if (details) cat( n," ", formatC(pwr, digits=6),"\n")
    if (iter>imax) break  
    # loop results in n with power too low
    # must step one up again. is done in the next loop
  }
  while (pwr<targetpower) {
    n    <- n+seqs
#    does not function stepsize too big in some cases    
#    if (abs(pwr-targetpower)>0.03) n  <- n+seqs
    iter <- iter+1
    fact <- bk/n
    sdm  <- sqrt(sd2*fact)
    df   <- eval(dfe)
    dfRR <- eval(dfRRe)
    
    p <- .power.scABEL(mlog, sdm, fact, sd2, df, sWR2, dfRR, nsims, 
                       ln_lBEL=log(theta1),ln_uBEL=log(theta2), 
                       CVswitch, r_const, CVcap, alpha=alpha)
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


