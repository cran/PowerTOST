# --------------------------------------------------------------------------
# power (or alpha) of 2-stage studies according to Potvin et. al. 
# methods "B" and "C", modified to include a futility criterion Nmax
# modified to use PE of stage 1 in sample size estimation
#
# author D.L.
# --------------------------------------------------------------------------
# require(PowerTOST)
# source("C:/Users/dlabes/workspace/PowerTOST/R/ss2x2.R")
# source("C:/Users/dlabes/workspace/PowerTOST/R/power.R")


power.2stage <- function(method=c("B","C"), alpha0=0.05, alpha=c(0.0294,0.0294),
                         n1, GMR, CV, pmethod=c("nct","exact"), targetpower=0.8, 
                         usePE=FALSE, theta0, theta1, theta2, Nmax=Inf, 
                         npct=c(0.05,0.25,0.5,0.75,0.95), nsims=1e5,
                         setseed=TRUE, print=TRUE, details=TRUE)
{
  if (missing(CV)) stop("CV must be given!")
  if (missing(n1)) stop("Number of subjects in stage 1 must be given!")
  
  if (missing(GMR)) GMR <- 0.95
  if (missing(theta1) & missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 <- 1/theta1
  
  if (missing(theta0)) theta0 <- GMR
  
  # check if Potvin B or C
  method  <- match.arg(method)
  # check if power calculation method is nct or exact
  pmethod <- match.arg(pmethod)
  
  ptm  <- proc.time()
  
  if (setseed) set.seed(1234567)

  ltheta1 <- log(theta1)
  ltheta2 <- log(theta2)
  lGMR    <- log(GMR)
  mlog    <- log(theta0)
  mse     <- CV2mse(CV)            
  BE      <- rep.int(NA, times=nsims)

  bk   <- 2   # 2x2x2 crossover design const
  
# ----- stage 1 ----------------------------------------------------------
  Cfact <- bk/n1
  df    <- n1-2
  tval  <- qt(1-alpha[1], df)
  sdm   <- sqrt(mse*Cfact)
  # simulate point est. via normal distribution
  pes   <- rnorm(n=nsims, mean=mlog, sd=sdm)
  # simulate mse via chi-squared distribution
  mses  <- mse*rchisq(n=nsims, df=df)/df
  # calculate power 
  if(method=="C"){
    # if method=C then calculate power for alpha0=0.05 and plan GMR
    pwr <- .calc.power(alpha=alpha0, ltheta1=ltheta1, ltheta2=ltheta2, 
                       diffm=lGMR, se=sqrt(mses), n=n1, df=df, bk=bk, 
                       method=pmethod)
    
    tval0 <- qt(1-alpha0, df)
    hw    <- tval0*sqrt(Cfact*mses)
    lower <- pes-hw
    upper <- pes+hw
    # fail or pass
    BE    <- lower>=ltheta1 & upper<=ltheta2
    # if power>0.8 then calculate CI for alpha=0.05
    # i.e. if power<0.8 then 
    BE[pwr<targetpower] <- NA # not yet decided
  }
  # method "B" or power<=0.8 in method "C"
  # calculate power for alpha=alpha[1]
  mses_tmp <- mses[is.na(BE)]
  pes_tmp  <- pes[is.na(BE)]
  BE1 <- rep.int(NA, times=length(mses_tmp))
  # calculate CI for alpha=alpha1
  hw    <- tval*sqrt(Cfact*mses_tmp)
  lower <- pes_tmp - hw
  upper <- pes_tmp + hw
  BE1   <- lower>=ltheta1 & upper<=ltheta2
  if (method=="C"){
    #if BE met -> PASS stop
    #if not BE -> goto sample size estimation i.e flag BE1 as NA
    BE1[!BE1] <- NA
  } else { 
    # method B
    # evaluate power at alpha[1]
    pwr <- .calc.power(alpha=alpha[1], ltheta1=ltheta1, ltheta2=ltheta2, 
                       diffm=lGMR, se=sqrt(mses_tmp), n=n1, df=df, bk=bk, 
                       method=pmethod)
    # if BE then decide BE regardless of power
    # if not BE and power<0.8 then goto stage 2
    BE1[ !BE1 & pwr<targetpower] <- NA 
  }
  # combine 'stage 0' from method C and stage 1
  BE[is.na(BE)] <- BE1
  # for calculating power in stage 1 below
  BE1 <- BE[!is.na(BE)]

  if(print & details){
    cat(nsims,"sims. Stage 1: Time consumed (min):\n")
    print(round((proc.time()-ptm)/60,2))
  }

  # ---------- stage 2 --------------------------------------------------
  ntot     <- rep(n1, times=nsims)
  # filter out those were stage 2 is necessary
  pes_tmp  <- pes[is.na(BE)]
  
  # Try to make all the stage 2 calculations without a loop
  # Maybe we are done with stage 1
  if(length(pes_tmp)>0){
    if(print & details){
      cat("Keep calm. Sample sizes for stage 2 will be estimated.\n")
      cat("May need some time. ")
    }
    mses_tmp <- mses[is.na(BE)]
    BE2      <- rep.int(NA, times=length(mses_tmp))
    # sample size for stage 2
    ptms <- proc.time()
    # use mse1 & pe1 if user decided so
    if (usePE){
      nt <- mapply(FUN=.sampleN.2x2, mse=mses_tmp, ltheta0=pes_tmp, 
                   MoreArgs=list(alpha=alpha[2], targetpower=targetpower, 
                                 ltheta1=ltheta1, ltheta2=ltheta2,
                                 method=pmethod))
    } else {
      nt <- mapply(FUN=.sampleN.2x2, mse=mses_tmp, 
                   MoreArgs=list(alpha=alpha[2], targetpower=targetpower, 
                                 ltheta0=lGMR, ltheta1=ltheta1, ltheta2=ltheta2,
                                 method=pmethod))
    }
    
    n2  <- ifelse(nt>n1, nt - n1, 0)
    # debug
    # cat("\n",length(nt[(n1+n2)==Nmax]), "values == Nmax\n")
    # print(table(as.factor(nt)))
    # debug print
    # print(summary(n2))
    
    #n2  <- ifelse(n2>100000, 100000, n2) # may not necessary
    
    if(print & details){
      cat("Time consumed (min):\n")
      print(round((proc.time()-ptms)/60,2))
    }
    # futility rule: if nt > Nmax -> stay with stage 1 result not BE
    # ntotal = n1 reasonable?
    if (is.finite(Nmax) | any(!is.finite(nt))){
      #sample size may return Inf if PE is used in ss estimation
      BE2[!is.finite(n2) | (n1+n2)>Nmax] <- FALSE
      # debug print
      # cat(sum(!BE2, na.rm=T)," cases with nt>Nmax or nt=Inf\n")
      # save the FALSE and NA in BE
      BE[is.na(BE)] <- BE2
      # filter out those were BE was yet not decided
      pes_tmp  <- pes_tmp[is.na(BE2)]
      mses_tmp <- mses_tmp[is.na(BE2)]
      n2       <- n2[is.na(BE2)]
    }
    m1    <- pes_tmp
    SS1   <- (n1-2)*mses_tmp
    nsim2 <- length(pes_tmp)
    # to avoid warnings for n2=0
    ow <- options("warn")
    options(warn=-1)
    m2    <- ifelse(n2>0, rnorm(n=nsim2, mean=mlog, sd=sqrt(mse*bk/n2)), 0)
    SS2   <- ifelse(n2>2, (n2-2)*mse*rchisq(n=nsim2, df=n2-2)/(n2-2), 0)
    # reset options
    options(ow) 
    SSmean <- ifelse(n2>0, (m1-m2)^2/(2/n1+2/n2), 0)
    nt    <- n1+n2
    df2   <- ifelse(n2>0, nt-3, n1-2)
    pe2   <- ifelse(n2>0, (n1*m1+n2*m2)/nt, pes_tmp)
    mse2  <- ifelse(n2>0, (SS1+SSmean+SS2)/df2, mses_tmp)
    # take care of memory
    rm(m1, m2, SS1, SS2, SSmean)
    # calculate CI for stage 2 with alpha[2]
    tval2 <- qt(1-alpha[2], df2)
    hw    <- tval2*sqrt(mse2*bk/nt)
    lower <- pe2 - hw
    upper <- pe2 + hw
    BE2   <- lower>=ltheta1 & upper<=ltheta2
    # combine stage 1 & stage 2
    ntot[is.na(BE)] <- nt
    BE[is.na(BE)]   <- BE2
    # debug print
    #cat("BE2", sum(BE2), "\n")
  } # end stage 2 calculations
  # output
  if (print) {
    if (details){
      cat("Total time consumed (min):\n")
      print(round((proc.time()-ptm)/60,2))
      cat("\n")
    }
    cat("Method ", method,":", sep="")
    if (method=="C") cat(" alpha0= ", alpha0, ",",sep="")
    cat(" alpha (s1/s2)=", alpha[1], alpha[2], "\n")
    cat("Futility criterion Nmax= ",Nmax,".\n", sep="")
    cat("CV= ",CV,"; n(stage 1)= ",n1,"; GMR= ",GMR, "\n", sep="")
    cat("BE margins = ",theta1," ... ", theta2,"\n", sep="")
    if (usePE) cat("Using PE and mse of stage 1 in sample size est.\n") else
      cat("Using GMR=",GMR, "and mse of stage 1 in sample size est.\n")
    cat("\n",nsims," sims at theta0= ", theta0, sep="")
    if (theta0<=theta1 | theta0>=theta2) cat(" (p(BE)='alpha').\n") else 
       cat(" (p(BE)='power').\n")
    cat("p(BE)   = ", sum(BE)/nsims,"\n", sep="")
    cat("p(BE) s1= ", sum(BE1)/nsims,"\n", sep="")
    cat("pct in stage 2= ", round(100*((nsims-length(BE1))/nsims),2),"%\n", sep="")
    cat("\nDistribution of N(total)\n")
    cat("- mean (range)= ", round(mean(ntot),1)," (", min(ntot)," ... ",
        max(ntot),")\n", sep="")
    cat("- percentiles\n")
    print(quantile(ntot, p=npct))
    cat("\n")
  } 
  #what shall we return?
  res <- list(method=method, alpha0=alpha0, alpha=alpha, CV=CV, n1=n1, GMR=GMR,
              theta1=theta1, theta2=theta2, theta0=exp(mlog), pmethod=pmethod,
              usePE=usePE, pBE=sum(BE)/nsims, pBEs1=sum(BE1)/nsims, 
              pct_s2=100*((nsims-length(BE1))/nsims), nmean=mean(ntot),
              nrange=range(ntot), nperc=quantile(ntot, p=npct), 
              nsims=nsims)
  if (print) return(invisible(res)) else return(res)
  
} #end function
