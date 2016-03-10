#-----------------------------------------------------------------------
# Iteratively adjust alpha to maintain the TIE <= nominal alpha for
# partial and full replicate design and scaled ABE via simulated
# (empirical) power
#
# Author: Helmut Schuetz
#-----------------------------------------------------------------------
scABEL.ad <-function(alpha = 0.05, theta0, theta1, theta2, CV = 0.3,
                     design = c("2x3x3", "2x2x4", "2x2x3"),
                     regulator = c("EMA", "ANVISA"), n, alpha.pre = 0.05,
                     imax=100, print = TRUE, details = FALSE, setseed = TRUE,
                     nsims = 1e6)
{
  ## Arguments:
  ##   alpha      Nominal alpha (in BE generally fixed to 0.05).
  ##              Lower value only if needed (e.g. to correct for
  ##              multiplicity).
  ##   theta0     If given, power is estimated for the expected GMR.
  ##              Note that the default is 0.90 (different from
  ##              sampleN.scABEL(), where the default is 0.95!
  ##   theta1     Lower margin. Defaults to 0.8.
  ##   theta2     Upper margin. Defaults to 1/theta1.
  ##   CV         Intra-subject CV(s) obtained in a replicate design.
  ##              (ratio, /not/ percent).
  ##              If given as a scalar, the CV of R.
  ##              If given as a vector, CV[1] /must/ be the CV of T and
  ##              CV[2] the CV of R. Important!
  ##              Defaults to 0.3 (maximum TIE for EMA).
  ##   design     "2x2x4", "2x2x3", "2x3x3". Defaults to "2x3x3".
  ##   regulator  "EMA" or "ANVISA". "ANVISA" requires extreme
  ##              adjustment close to CVwR 40%. Realistic?
  ##              Cave: ANVISA's requirements are unofficial.
  ##   n          Total sample size or a vector of subjects/sequences.
  ##   nsims      Simulations for the TIE. Should not be <1e6.
  ##   imax       max. steps in sample size search
  ##   print      Boolean (FALSE returns a list of results).
  ##   details    Boolean (runtime, number of simulations).
  ##   alpha.pre  Pre-specified level.
  ##   setseed    Boolean (default TRUE uses set.seed(123456)).
  ## Returns:
  ##   alpha.adj  Iteratively adjusted alpha which does not inflate the
  ##              TIE (for given CVwR and n).
  ##   CI.adj     The adjusted confidence interval in percent, where
  ##              CI.adj = 100(1-2*alpha.adj).
  ##   TIE.unadj  The empiric Type I Error based on nominal alpha.
  ##   TIE.adj    TIE based on adjusted alpha.
  ##   rel.change Relative change in risk (%) compared to nominal alpha.
  ##   If theta0 is given:
  ##   pwr.unadj  Power for alpha (or, if given, alpha.pre).
  ##   pwr.adj    Power for adjusted alpha.
  ##   rel.loss   Relative loss in power if the sample size was planned
  ##              for alpha and will be evaluated with alpha.adj,
  ##              where rel.loss = 100(pwr.adj ? pwr.unadj)/pwr.unadj
  ##   If alpha.pre is given:
  ##   Assessment of TIE; alpha.pre is justified if n.s. > alpha.
  ######################################################################
  ## Tested on Win 7 Pro SP1 64bit
  ##   R 3.2.3 64bit (2015-12-10), PowerTOST 1.3-02 (2015-12-02)
  ######################################################################
  env <- as.character(Sys.info()[1]) # get info about the OS
  if ((env == "Windows") || (env == "Darwin")) flushable <- TRUE
    else flushable <- FALSE # supress flushing on other OS's
  # acceptance range defaults
  if (missing(theta1) && missing(theta2)) theta1 <- 0.8
  if (missing(theta2)) theta2 = 1/theta1
  # check theta0
  if (missing(theta0)) theta0 <- 0.9
  if (theta0 < theta1 || theta0 > theta2)
    stop("theta0 must be within [theta1, theta2]")
  # check regulator arg
  regulator <- toupper(regulator)
  regulator <- match.arg(regulator)
  # set iteration tolerance for uniroot(). Must be higher for ANVISA.
  if (regulator == "EMA") tol <- 1e-5 else tol <- 1e-6
  design <- match.arg(design)
  CVwT <- CV[1]
  if (length(CV) == 2) CVwR <- CV[2] else CVwR <- CVwT
  no <- 0 # simulation counter
  if (details) ptm <- proc.time()
  # Was geht hier ab? Sollten wir unbedingt in der man page beschreiben.
  if (missing(n) || is.na(n)) {
    if (is.na(alpha.pre) || (alpha.pre != alpha)) {
      al <- alpha.pre    # If pre-specified, use alpha.pre
    } else {
      al <- alpha        # If not, use alpha (commonly 0.05)
    }
    n <- sampleN.scABEL(alpha = al, CV = CV, theta0 = theta0,
                        design = design, regulator = regulator,
                        imax=imax, print = FALSE, details = FALSE,
                        setseed = setseed)[["Sample size"]]
    if (is.na(n))
      stop(paste0("Sample size search in sampleN.scABEL() failed.",
                  "\nRestart with an explicit high n (>1000)."))
    if (sum(n) < 6) stop("Sample size too low.")
    no <- 1e5
  }
  if (alpha.pre > alpha) {
    warning(paste0("alpha.pre > alpha doesn't make sense.",
                   "\nalpha.pre was set to alpha."))
    alpha.pre <- alpha
  }
  seqs <- as.numeric(substr(design, 3, 3)) # subjects / sequence
  if (length(n) == 1) n <- nvec(n, seqs) # vectorize n

  # los gehts
  pwr       <- rep(NA, 2) # initialize vectors: pwr
  TIE       <- rep(NA, 2) # TIE
  alpha.adj <- NA   # adjusted alpha
  opt <- function(x) power.scABEL(alpha = x, CV = CV, theta0 = U, n = n,
                                  regulator = regulator, design = design,
                                  nsims = nsims, setseed = setseed) - alpha
  # Finds adjusted alpha which gives TIE as close as possible to alpha.
  sig  <- binom.test(x = round(alpha*nsims, 0), n = nsims,
                     alternative = "less", conf.level = 1 - alpha)$conf.int[2]
  method <- "ABE"
  if ((regulator == "EMA" && CVwR > 0.3) ||
      (regulator == "ANVISA" && CVwR > 0.4)) method <- "ABEL"
  U <- scABEL(CV = CVwR, regulator = regulator)[["upper"]]
  # Simulate at the upper (expanded) limit. For CVwR 30% that's 1.25.
  # Due to the symmetry simulations at the lower limit (0.8) would work
  # as well.
  if (alpha.pre != alpha) {
    al <- alpha.pre # If pre-specified, use alpha.pre.
  } else {
    al <- alpha     # If not, use alpha (commonly 0.05).
  }
  designs <- c("2x2x4", "2x2x3", "2x3x3")
  type    <- c("RTRT|TRTR", "RTR|TRT", "RRT|RTR|TRR") # clear words
  if (print) { # Show input to keep the spirits of the user high.
    cat("\n+++++++++++ scaled (widened) ABEL +++++++++++\n")
    cat("         iteratively adjusted alpha\n")
    cat("---------------------------------------------\n")
    cat("Study design: ")
    cat(paste0(design, " (", type[match(design, designs)], ")\n"))
    cat("log-transformed data (multiplicative model)\n")
    cat(formatC(nsims, format = "d", big.mark = ",", decimal.mark = "."),
        "studies in each iteration simulated.\n\n")
    txt <- paste0("CVwR ", sprintf("%.4g", CVwR))
    if (length(CV) == 2) {
      txt <- paste0(txt, ", CVwT ", sprintf("%.4g", CVwT), ", ")
    } else {
      txt <- paste0(txt, ", ")
    }
    cat(paste0(txt, "n(i) ", paste0(n, collapse = "|"), " (N ", sum(n),
                    ")\n"))
    txt <- paste0("Nominal alpha                 : ", signif(alpha, 5))
    if (!is.na(alpha.pre) && (alpha.pre != alpha)) {
      txt <- paste0(txt, ", pre-specified alpha ", alpha.pre, "\n")
    } else {
      txt <- paste(txt, "\n")
    }
    cat(txt)
    cat("Null (true) ratio             :", sprintf("%.3f", theta0), "\n")
    cat(paste0("Regulatory settings           : ", regulator, " (",
               method, ")\n"))
    cat(paste0("Significance limit of TIE     : ", signif(sig, 5), "\n"))
    if (flushable) flush.console() # advance console output.
  }
  TIE[1] <- power.scABEL(alpha = al, CV = CV, theta0 = U, n = n,
                         design = design, regulator = regulator,
                         nsims = nsims, setseed = setseed)
  no <- no + nsims
  pwr[1] <- power.scABEL(alpha = al, CV = CV, theta0 = theta0,
                         n = n, design = design, regulator = regulator,
                         setseed = setseed)
  no <- no + 1e5
  if (TIE[1] > sig) { # adjust only if needed (significant inflation)
    x         <- uniroot(opt, interval = c(0, alpha), tol = tol)
    alpha.adj <- x$root
    TIE[2] <- power.scABEL(alpha = alpha.adj, CV = CV, theta0 = U, n = n,
                           design = design, regulator = regulator,
                           nsims = nsims, setseed = setseed)
    pwr[2] <- power.scABEL(alpha = alpha.adj, CV = CV, theta0 = theta0,
                           n = n, design = design, regulator = regulator,
                           setseed = setseed)
  }
  if (!is.na(alpha.adj)) no <- no + nsims*x$iter
  if (details) run.time <- proc.time() - ptm
  if (print) { # fetch and print results
    txt <- paste0("Empiric TIE for alpha ", sprintf("%.4f", al), "  : ",
                  sprintf("%.4f", TIE[1]))
    if (TIE[1] > sig || alpha.pre != alpha) {
      rel.change <- 100*(TIE[1] - alpha)/alpha
      if (details) {
        txt <- paste0(txt, " (rel. change of risk: ",
                      sprintf("%+1.3g%%", rel.change), ")")
      }
    }
    if (!is.na(pwr[1])) {
      pwr.unadj <- pwr[1]
      txt <- paste0(txt, "\nPower for theta0 ", sprintf("%.3f", theta0),
                         "        : ", sprintf("%.3f", pwr.unadj))
    }
    if (TIE[1] > sig) {
      if (alpha.adj >= 0.01) {
        txt <- paste0(txt, "\nIteratively adjusted alpha    : ",
                      sprintf("%.4f", alpha.adj),
                      "\nEmpiric TIE for adjusted alpha: ",
                      sprintf("%.4f", TIE[2]))
      } else {
        txt <- paste0(txt, "\nIteratively adjusted alpha    : ",
                      signif(alpha.adj, 3),
                      "\nEmpiric TIE for adjusted alpha: ",
                      sprintf("%.4f", TIE[2]))
      }
    if (!is.na(pwr[2])) {
      pwr.adj <- pwr[2]
      txt <- paste0(txt, "\nPower for theta0 ",sprintf("%.3f", theta0),
                         "        : ", sprintf("%.3f", pwr.adj))
      if (details) {
        txt <- paste0(txt, " (rel. impact: ", sprintf("%+1.3g%%",
                      100*(pwr[2] - pwr[1])/pwr[1]), ")")
      }
      txt <- paste(txt, "\n\n")
    } else {
      txt <- paste0(txt, "\n\n")
    }
      if (details) {
        txt <- paste0(txt, "Runtime    : ", signif(run.time[3], 3),
                      " seconds\nSimulations: ",
                      formatC(no, format = "d", big.mark = ",",
                        decimal.mark = "."), " (",
                      (no - no %% nsims - nsims) / nsims,
                      " iterations)\n\n")
      }
    } else {
      txt <- paste0(txt, "\nNo significant inflation of the TIE ",
                    "expected; ")
      ifelse(alpha.pre == alpha,
        txt <- paste0(txt, "no adjustment of alpha is required.\n\n"),
        txt <- paste0(txt, "the chosen pre-specified alpha is ",
                      "justified.\n\n"))
    }
  cat(txt)
  } else { # Prepare and return list of results.
    res <- list(regulator = regulator, method = "ABEL", design = design,
                type = type[match(design, designs)], alpha = alpha,
                alpha.pre = alpha.pre, CV = CV, n = sum(n), theta0 = theta0,
                TIE.unadj = signif(TIE[1], 5),
                rel.change = ifelse(!is.na(alpha.adj),
                                  signif(100*(TIE[1] - alpha)/alpha, 5), NA),
                pwr.unadj = signif(pwr[1], 5), alpha.adj = signif(alpha.adj, 5),
                TIE.adj = signif(TIE[2], 5), pwr.adj = signif(pwr[2], 5),
                rel.loss = ifelse(!is.na(pwr[2]),
                                signif(100*(pwr[2] - pwr[1])/pwr[1], 5), NA),
                sims = no)
    return(res)
  }
}
# Examples
# 1. using all defaults
#   scABEL.ad()
# should return:
#   +++++++++++ scaled (widened) ABEL +++++++++++
#            iteratively adjusted alpha
#   ---------------------------------------------
#   Study design: 2x3x3 (RRT|RTR|TRR)
#   log-transformed data (multiplicative model)
#   1,000,000 studies in each iteration simulated.
#
#   CVwR 0.3, n(i) 18|18|18 (N 54)
#   Nominal alpha                 : 0.05
#   Null (true) ratio             : 0.900
#   Regulatory settings           : EMA (ABE)
#   Significance limit of TIE     : 0.05036
#   Empiric TIE for alpha 0.0500  : 0.0719
#   Iteratively adjusted alpha    : 0.0339
#   Empiric TIE for adjusted alpha: 0.0500
#   Power for theta0 0.900        : 0.764
#
# 2. Explore the impact on power.
#   scABEL.ad(regulator="EMA", design="2x2x4", CV=0.3, n=34, details=TRUE, theta0=0.9)
# should return:
#   +++++++++++ scaled (widened) ABEL +++++++++++
#            iteratively adjusted alpha
#   ---------------------------------------------
#   Study design: 2x2x4 (RTRT|TRTR)
#   log-transformed data (multiplicative model)
#   1,000,000 studies in each iteration simulated.
#
#   CVwR 0.3, n(i) 17|17 (N 34)
#   Nominal alpha                 : 0.05
#   Null (true) ratio             : 0.900
#   Regulatory settings           : EMA (ABE)
#   Significance limit of TIE     : 0.05036
#   Empiric TIE for alpha 0.0500  : 0.0816 (rel. change of risk: +63.3%)
#   Power for theta0 0.900        : 0.803
#   Iteratively adjusted alpha    : 0.0286
#   Empiric TIE for adjusted alpha: 0.0500
#   Power for theta0 0.900        : 0.725 (rel. impact: -9.68%)
#
#   Runtime    : 5.46 seconds
#   Simulations: 5,100,000 (4 iterations)
#
# 3. Explore whether a pre-specified alpha maintains the consumer's risk;
#    different CVs of T and R and slight unbalance in the data expected.
#   scABEL.ad(regulator="EMA", design="2x2x4", CV=CVp2CV(0.3, ratio=0.75), n=c(15, 14), alpha.pre=0.025)
# should return:
#   +++++++++++ scaled (widened) ABEL +++++++++++
#            iteratively adjusted alpha
#   ---------------------------------------------
#   Study design: 2x2x4 (RTRT|TRTR)
#   log-transformed data (multiplicative model)
#   1,000,000 studies in each iteration simulated.
#
#   CVwR 32.17%, CVwT 27.69%, n(i) 15|14 (N 29)
#   Nominal alpha                 : 0.05, pre-specified alpha 0.025
#   Null (true) ratio             : 0.900 (default)
#   Regulatory settings           : EMA (ABE)
#   Significance limit of TIE     : 0.05036
#   Empiric TIE for alpha 0.0250  : 0.0394
#   No significant inflation of the TIE expected; the chosen pre-specified alpha is justified.#
#
# 4. Assign to a variable and subsequently call the results
#   x <- scABEL.ad(regulator="EMA", design="2x2x4", CV=0.3, n=34, theta0=0.9, print=FALSE, alpha.pre=0.025)
# TIE for the requested alpha
#   x$TIE.unadj
#   [1] 0.044427 # No inflation of the TIE; 0.025 is justified.
#   x$pwr.unadj
#   [1] 0.70537  # Consider to increase the sample size!
