#---------------------------------------------------------------------------
# Simulate partial and full replicate design and scaled ABE power
# 
# Author: dlabes
#---------------------------------------------------------------------------

# degrees of freedom for the RR  analysis: 
# model with subj, period like EMA Q&A
# the inclusion of sequence changes only the distribution of df's
# between seq df=seq-1 and sub(seq) df=n-1 - seq-1. 
# 2x3x3  2*n measurements df = 2*n-1
#        n subjects       df = n-1
#        3 periods        df = 2
#                   ->  dfRR = n-2
# 2x2x4  2*n measurements df = 2*n-1
#        n subjects       df = n-1
#        4 periods        df = 3
#                   ->  dfRR = n-3
# But cave! The EMA set I has only 2 df for period. Due to imbalance?
# Another possibility is using the contrasts R-R and analyze by sequence. 
# Then the dfRR = n-seq.
# 2x3x3  dfRR = n-3
# 2x2x4  dfRR = n-2


power.scABEL <- function(alpha=0.05, theta1, theta2, theta0, CV, n,   
                         design=c("2x3x3", "2x2x4"), regulator=c("EMA","FDA"),
                         nsims=1E5, details=FALSE, setseed=TRUE)
{
  if (missing(CV)) stop("CV must be given!")
  if (missing(n))  stop("Number of subjects n must be given!")

  if (missing(theta0)) theta0 <- 0.95
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 <- 1/theta1
  
  ptm <- proc.time()
  
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
  
  design <- match.arg(design)
  if (design=="2x3x3") {
    seqs <- 3
    bkni <- 1/6
    bk   <- 1.5
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
    seqs <- 2
    bkni <- 1/4
    bk   <- 1
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
  CVwT <- CV[1]
  if (length(CV)==2) CVwR <- CV[2] else CVwR <- CVwT
  s2WT <- log(1.0 + CVwT^2)
  s2WR <- log(1.0 + CVwR^2)
  # sd^2 (variance) of the differences T-R from their components
  sd2  <- (sD2 + (s2WT + s2WR)/2) # is this correct for partial replicate?

  # sd^2 of the differences T-R from their components
  sd2  <- (sD2 + (s2WT + s2WR)/2) # is this correct for partial replicate?
  
  if (length(n)==1){
    # for unbalanced designs we divide the ns by ourself
    # to have only little imbalance
    ni <- round(n/seqs,0)
    nv <- rep.int(ni, times=seqs-1)
    nv <- c(nv, n-sum(nv))
    if (nv[length(nv)]!=ni){
      message("Unbalanced design. n(i)=", paste(nv, collapse="/"),
              " assumed.")
    } 
    fact <- sum(1/nv)*bkni
    n <- sum(nv)
  } else {
    # check length
    if (length(n)!=seqs) stop("n must be a vector of length=",seqs,"!")
    
    fact <- sum(1/n)*bkni
    n <- sum(n)
  }
  # sd of the sample mean T-R (point estimator)
  sdm  <- sqrt(sd2*fact)
  mlog <- log(theta0)
  df   <- eval(dfe)
  dfRR <- eval(dfRRe)
  
  if(setseed) set.seed(123456)
  p <- .power.scABEL(mlog, sdm, fact, sd2, df, s2WR, dfRR, nsims, 
                     ln_lBEL=log(theta1),ln_uBEL=log(theta2), 
                     CVswitch, r_const, CVcap, alpha=alpha)
    
  if (details) {
    cat(nsims,"sims. Time elapsed (sec):\n")
    print(proc.time()-ptm)
    cat("p(BE-ABE)=", p["BEabe"],"; p(BE-wABEL)=", p["BEwl"],
        "; p(BE-PE)=", p["BEpe"],"\n\n")
  }
  # return the 'power'
  as.numeric(p["BE"])
}

.power.scABEL <- function(mlog, sdm, fact, sd2, df, s2WR, dfRR, nsims, 
                           CVswitch=0.3, r_const=0.760, CVcap=0.5,
                           ln_lBEL=log(0.8), ln_uBEL=log(1.25), alpha=0.05)
{
  tval <- qt(1-alpha,df)
  
  counts <- rep.int(0, times=4)
  names(counts) <- c("BE", "BEpe", "BEwl","BEabe")
  # to avoid memory problems
  chunks <- 1
  nsi    <- nsims
  if (nsims>1E7) {
    chunks <- round(nsims/1E7,0)
    nsi    <- 1E7
  } 
  for (iter in 1:chunks){
    # simulate sample mean via its normal distribution
    means  <- rnorm(nsi, mean=mlog, sd=sdm)
    # simulate sample sd2 via chi-square distri
    sd2s   <- sd2*rchisq(nsi, df)/df
    # simulate sample value s2WR via chi-square distri
    s2WRs  <- s2WR*rchisq(nsi, dfRR)/dfRR
    CVwr   <- sqrt(exp(s2WRs)-1)
    # EMA limits in log-domain
    lABEL   <- -sqrt(s2WRs)*r_const
    uABEL   <- +sqrt(s2WRs)*r_const
    lABEL[CVwr<=CVswitch] <- ln_lBEL
    uABEL[CVwr<=CVswitch] <- ln_uBEL
    # cap
    if (is.finite(CVcap)){
      lABEL[CVwr>0.5] <- -sqrt(log(1.0 + CVcap^2))*r_const
      uABEL[CVwr>0.5] <- +sqrt(log(1.0 + CVcap^2))*r_const
    }
    # 90% CIs for T-R
    hw  <- tval*sqrt(sd2s*fact)
    lCL <- means - hw 
    uCL <- means + hw
    rm(hw)
    # conventional ABE
    BEABE <- ( ln_lBEL<=lCL & uCL<=ln_uBEL )
    # 90% CI in widened limits? 
    BE   <- (lABEL<=lCL & uCL<=uABEL)
    # point est. constraint true?
    BEpe <- ( means>=ln_lBEL & means<=ln_uBEL )
    
    counts["BEabe"] <- counts["BEabe"] + sum(BEABE)
    counts["BEpe"]  <- counts["BEpe"]  + sum(BEpe)
    counts["BEwl"]  <- counts["BEwl"]  + sum(BE)
    counts["BE"]    <- counts["BE"]    + sum(BE & BEpe)
  } # end over chunks
  # return the counts
  counts/nsims
}