#---------------------------------------------------------------------------
# Simulate full replicate design and calculate scaled ABE power
# according to FDA Warfarin guidance
#
# Author: dlabes
#---------------------------------------------------------------------------
# 2x2x4  dfRR = n-2

power.NTIDFDA <- function(alpha=0.05, theta1, theta2, theta0, CV, n,   
                          nsims=1E5, details=FALSE, setseed=TRUE)
{
  if (missing(CV)) stop("CV must be given!", call.=FALSE)
  if (missing(n))  stop("Number of subjects n must be given!", call.=FALSE)
   
  if (missing(theta0)) theta0 <- 0.975  # tighter content limits for NTID
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 <- 1/theta1
  
  ptm <- proc.time()
  
  # for later enhancement taking into account the 
  # subject-by-formulation interaction
  s2D  <- 0 
  CVwT <- CV[1]
  if (length(CV)==2) CVwR <- CV[2] else CVwR <- CVwT
  s2wT <- log(1.0 + CVwT^2)
  s2wR <- log(1.0 + CVwR^2)

  # FDA constant
  r_const  <- -log(0.9)/0.10
  
  # design only 2x2x4
  seqs  <- 2
  dfe   <- parse(text="n-2", srcfile=NULL)
  dfRRe <- parse(text="n-2", srcfile=NULL)
  # sd^2 of the differences T-R from their components
  Emse  <- (s2D + (s2wT + s2wR)/2) 
  
  if (length(n)==1){
    # we assume n=ntotal
    # for unbalanced designs we divide the ns by ourself
    # in such a way that we have only small imbalance
    nv <- nvec(n=n, grps=seqs)
    if (nv[1]!=nv[length(nv)]){
      message("Unbalanced design. n(i)=", paste(nv, collapse="/"), " assumed.")
    }
    C3 <- sum(1/nv)/seqs^2
    n  <- sum(nv)
  } else {
    # we assume n = vector of n's in sequences
    # check length
    if (length(n)!=seqs) stop("n must be a vector of length=",seqs,"!", call.=FALSE)
    C3 <- sum(1/n)/seqs^2
    n  <- sum(n)
  }
  # sd of the mean T-R (point estimator)
  sdm  <- sqrt(Emse*C3)
  mlog <- log(theta0)
  df   <- eval(dfe)
  dfRR <- eval(dfRRe)
  
  dfTT <- dfRR       # at least for the 2x2x4 design
 
  if(setseed) set.seed(123456)
  
  p <- .power.NTID(mlog, sdm, C3, Emse, df, s2wR, dfRR, s2wT, dfTT, nsims, 
                   r_const, ln_lBEL=log(theta1),ln_uBEL=log(theta2), alpha=alpha)
    
  if (details) {
    cat(nsims,"sims. Time elapsed (sec):\n")
    print(proc.time()-ptm)
    cat("p(BE-ABE)=", p["BEabe"],"; p(BE-SABEc)=", p["BEul"],
        "; p(BE-sratio)=", p["BEsratio"],"\n")
    cat("\n")
  }
  # return the 'power'
  as.numeric(p["BE"])
}

# working horse of RSABE for NTID's
.power.NTID <- function(mlog, sdm, C3, Emse, df, s2wR, dfRR, s2wT, dfTT, nsims, 
                        r_const=-log(0.9)/0.1, ln_lBEL=log(0.8), ln_uBEL=log(1.25), 
                        alpha=0.05)
{
  tval     <- qt(1-alpha,df)
  chisqval <- qchisq(1-alpha, dfRR)
  r2const  <- r_const^2
  Fval     <- qf(1-alpha, dfTT, dfRR, lower.tail=FALSE)
  
  counts   <- rep.int(0, times=4)
  names(counts) <- c("BE", "BEsratio", "BEul","BEabe")
  # to avoid memory problems for high number of sims
  # we are working with chunks of 1e7
  chunks   <- 1
  nsi      <- nsims
  if (nsims>1E7) {
    chunks <- round(nsims/1E7,0)
    nsi    <- 1E7
  } 
  for (iter in 1:chunks) {
    # simulate sample mean via its normal distribution
    means  <- rnorm(nsi, mean=mlog, sd=sdm)
    # simulate sample sd2s via chi-square distri
    sd2s   <- Emse*C3*rchisq(nsi, df)/df
    # simulate sample value s2wRs via chi-square distri
    s2wRs  <- s2wR*rchisq(nsi, dfRR)/dfRR
    # simulate sample value s2wTs via chi-square distri
    s2wTs  <- s2wT*rchisq(nsi, dfTT)/dfTT
    
    SEs <- sqrt(sd2s)
    
    # conventional (1-2*alpha) CI's for T-R
    hw  <- tval*SEs
    lCL <- means - hw 
    uCL <- means + hw
    # conventional ABE
    BEABE     <- ((ln_lBEL<=lCL) & (uCL<=ln_uBEL))
    
    # upper 95% CI linearized SABE criterion
    # with -SEs^2 the 'unknown' x from the warfarin guidance
    Em <- means^2 - SEs^2  
    Es <- r2const*s2wRs
    #Cm <- (abs(means) + hw)^2
    Cm <- ifelse(abs(lCL)>abs(uCL),abs(lCL)^2,abs(uCL)^2)
    Cs <- Es*dfRR/chisqval    
    SABEc95 <- Em - Es + sqrt((Cm-Em)^2 + (Cs-Es)^2)
    BEscABE   <- (SABEc95 <= 0)
    # save memory
    rm(SEs, hw, Em, Es, Cm, Cs)
    
    # upper limit of ratio swT/swR
    ul_sratio <- sqrt(s2wTs/s2wRs/Fval)
    # upper limit <= 2.5?
    BEsratio  <- ul_sratio <= 2.5
    
    # debug print
    if (nsims<=50){
      print(head(data.frame(pe=means,lCL=lCL,uCL=uCL, BEABE, s2wT=s2wTs, s2wR=s2wRs, 
                            SABEc95, BEscABE, sratio=sqrt(s2wTs/s2wRs), ul_sratio, 
                            BEsratio), n=50))
    } 
    
    counts["BEabe"]    <- counts["BEabe"]    + sum(BEABE)
    counts["BEul"]     <- counts["BEul"]     + sum(BEscABE)
    counts["BEsratio"] <- counts["BEsratio"] + sum(BEsratio)
    # test without s-ratio
    #counts["BE"]       <- counts["BE"]       + sum(BEscABE & BEABE)
    counts["BE"]       <- counts["BE"]       + sum(BEscABE & BEABE & BEsratio)
    
  } # end over chunks
  # return the pBEs
  counts/nsims
}