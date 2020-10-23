## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#"
)

## ----setup--------------------------------------------------------------------
library(PowerTOST) # attach the library

## ----example1a----------------------------------------------------------------
sampleN.TOST(CV = 0.30)

## ----example1b----------------------------------------------------------------
sampleN.TOST(CV = 0.30, details = FALSE, print = FALSE)[["Sample size"]]

## ----example1c----------------------------------------------------------------
power.TOST(CV = 0.30, n = 39)

## ----example1d----------------------------------------------------------------
designs <- c("2x2x2", "2x2x3", "2x3x3", "2x2x4")
# data.frame of results
res <- data.frame(design = designs, n = NA, power = NA, n.do = NA,
                  power.do = NA, stringsAsFactors = FALSE)
for (i in 1:4) {
  # print = FALSE and details = FALSE suppress output to the console
  # we are only interested in columns 7-8
  # let's also calculate power for one dropout
  res[i, 2:3] <- sampleN.TOST(CV = 0.30, design = res$design[i],
                              print = FALSE, details = FALSE)[7:8]
  res[i, 4]   <-  res[i, 2] - 1
  res[i, 5]   <- suppressMessages(
                    power.TOST(CV = 0.30, design = res$design[i],
                               n = res[i, 4]))
}
print(res, row.names = FALSE)

## ----example1e, echo = FALSE--------------------------------------------------
designs <- c("2x2x2", "2x2x3", "2x3x3", "2x2x4")
res <- data.frame(design = rep(NA, 4), name = NA, n = NA, formula = NA,
                  df = NA, t.value = NA, stringsAsFactors = FALSE)
res[, c(2:1, 4)] <- known.designs()[which(known.designs()[, 2] %in% designs),
                                    c(9, 2:3)]
for (i in 1:4) {
  res$n[i]  <- sampleN.TOST(CV = 0.30, design = res$design[i],
                            print = FALSE, details = FALSE)[["Sample size"]]
  e         <- parse(text=res[i, 4], srcfile=NULL)
  n         <- res$n[i]
  res[i, 5] <- eval(e)
  res$t.value[i] <- signif(qt(1-0.05, df = res[i, 5]), 4)
}
res <- res[with(res, order(-n, design)), ]
print(res, row.names = FALSE)

## ----example1f----------------------------------------------------------------
grouping <- function(capacity, n) {
  # split sample size into >=2 groups based on capacity
  if (n <= capacity) { # make equal groups
    ngrp <- rep(ceiling(n/2), 2)
  } else {             # at least one = capacity
    ngrp    <- rep(0, ceiling(n / capacity))
    grps    <- length(ngrp)
    ngrp[1] <- capacity
    for (j in 2:grps) {
      n.tot <- sum(ngrp) # what we have so far
      if (n.tot + capacity <= n) {
        ngrp[j] <- capacity
      } else {
        ngrp[j] <- n - n.tot
      }
    }
  }
  return(ngrp = list(grps = length(ngrp), ngrp = ngrp))
}
CV        <- 0.30
capacity  <- 24 # clinical capacity
res       <- data.frame(n = NA, grps = NA, n.grp = NA, m.1 = NA, m.2 = NA)
x         <- sampleN.TOST(CV = CV, print = FALSE, details = FALSE)
res$n     <- x[["Sample size"]]
res$m.1   <- x[["Achieved power"]]
x         <- grouping(capacity = capacity, n = res$n)
res$grps  <- x[["grps"]]
ngrp      <- x[["ngrp"]]
res$n.grp <- paste(ngrp, collapse = "|")
res$m.2   <- power.TOST.sds(CV = CV, n = res$n, grps = res$grps,
                            ngrp = ngrp, gmodel = 2, progress = FALSE)
res$loss <- 100*(res$m.2 - res$m.1)/res$m.1
names(res)[2:6] <- c("groups", "n/group", "pooled model",
                     "group model", "loss (%)")
res[1, 4:6] <- sprintf("%6.4f", res[1, 4:6])
cat("Achieved power and relative loss\n"); print(res, row.names = FALSE)

## ----example2-----------------------------------------------------------------
sampleN.RatioF(CV = 0.20, CVb = 0.40)

## ----example3a----------------------------------------------------------------
sampleN.TOST(CV = 0.20, theta0 = 0.92)

## ----example3b----------------------------------------------------------------
df <- 16 - 2 # degrees of freedom of the 2x2x2 crossover pilot
CVCL(CV = 0.20, df = df, side = "upper", alpha = 0.20)[["upper CL"]]

## ----example3c----------------------------------------------------------------
CL.upper <- CVCL(CV = 0.20, df = 16 - 2, side = "upper",
                 alpha = 0.20)[["upper CL"]]
res <- sampleN.TOST(CV = CL.upper, theta0 = 0.92, print = FALSE)
print(res[7:8], row.names = FALSE)

## ----example3d----------------------------------------------------------------
CL.upper <- CVCL(CV = 0.20, df = 16 - 2, side = "upper",
                 alpha = 0.20)[["upper CL"]]
power.TOST(CV = CL.upper, theta0 = 0.92, n = 28)

## ----example3e----------------------------------------------------------------
power.TOST(CV = 0.22, theta0 = 0.90, n = 40)

## ----example3f----------------------------------------------------------------
expsampleN.TOST(CV = 0.20, theta0 = 0.92, prior.type = "CV",
                prior.parm = list(m = 16, design = "2x2x2"))

## ----example3g----------------------------------------------------------------
expsampleN.TOST(CV = 0.20, theta0 = 0.92, prior.type = "theta0",
                prior.parm = list(m = 16, design = "2x2x2"))

## ----example3h----------------------------------------------------------------
expsampleN.TOST(CV = 0.20, theta0 = 0.92, prior.type = "both",
                prior.parm = list(m = 16, design = "2x2x2"),
                details = FALSE)

## ----example4-----------------------------------------------------------------
CV  <- 0.214
res <- data.frame(target = c(rep(0.8, 5), rep(0.9, 5)),
                  theta0 = rep(c(1, seq(0.95, 0.92, -0.01)), 2),
                  n.1 = NA, power = NA,
                  sigma.u = rep(c(0.001, seq(0.05, 0.08, 0.01)), 2),
                  n.2 = NA, assurance = NA)
for (i in 1:nrow(res)) {
  res[i, 3:4] <- sampleN.TOST(CV = CV, targetpower = res$target[i],
                              theta0 = res$theta0[i], details = FALSE,
                              print = FALSE)[7:8]
  res[i, 6:7] <- expsampleN.TOST(CV = CV, targetpower = res$target[i],
                                 theta0 = 1, prior.type = "theta0",
                                 prior.parm = list(sem = res$sigma.u[i]),
                                 details = FALSE, print = FALSE)[9:10]
}
res                 <- signif(res, 3)
res[, 5]            <- sprintf("%.2f", res[, 5])
names(res)[c(3, 6)] <- "n"
print(res, row.names = FALSE)

## ----example4a----------------------------------------------------------------
planned  <- "2x2x2"
logscale <- FALSE
theta0   <- -5
theta1   <- -15
theta2   <- +15 # if not given, -theta1 is used
SD.resid <-  20 # residual standard deviation
sampleN.TOST(CV = SD.resid, theta0 = theta0, theta1 = theta1,
             theta2 = theta2, logscale = logscale, design = planned)

## ----example4b----------------------------------------------------------------
known    <- known.designs()[, c(2, 6)]           # extract relevant information
bk       <- known[known$design == planned, "bk"] # retrieve design constant
cat(paste0("The design constant for design \"", planned, "\" is ", bk, "\n"))
SD.delta <-  28                                  # standard deviation of the difference
SD.resid <-  SD.delta / sqrt(bk)                 # convert to residual SD
sampleN.TOST(CV = SD.resid, theta0 = theta0, theta1 = theta1,
             theta2 = theta2, logscale = logscale, design = planned)

## ----higher_order1------------------------------------------------------------
CV  <- 0.20
res <- data.frame(design = c("3x6x3", "2x2x2"), n = NA, power = NA,
                  stringsAsFactors = FALSE)
for (i in 1:2) {
  res[i, 2:3] <- sampleN.TOST(CV = CV, design = res$design[i],
                              details = FALSE, print = FALSE)[7:8]
}
print(res, row.names = FALSE)

## ----higher_order2------------------------------------------------------------
CV  <- 0.20
res <- data.frame(design = c("4x4", "2x2x2"), n = NA, power = NA,
                  stringsAsFactors = FALSE)
for (i in 1:2) {
  res[i, 2:3] <- sampleN.TOST(CV = CV, design = res$design[i],
                              details = FALSE, print = FALSE)[7:8]
}
print(res, row.names = FALSE)

## ----power1-------------------------------------------------------------------
n <- c(14 - 1, 14 - 2) # 14 dosed in each sequence
round(100*CI.BE(pe = 0.90, CV = 0.25, n = n), 2)

## ----power2-------------------------------------------------------------------
power.TOST(CV = 0.25, theta0 = 0.90, n = c(13, 12)) # observed values

## ----power3-------------------------------------------------------------------
power.TOST(CV = 0.20, theta0 = 0.92, n = 28)     # assumed in planning
power.TOST(CV = 0.25, theta0 = 1, n = c(13, 12)) # observed but wrong T/R-ratio

## ----pooling------------------------------------------------------------------
CVs <- ("  CV |  n | design | study
         0.20 | 16 |  2x2x2 | pilot
         0.25 | 25 |  2x2x2 | pivotal")
txtcon <- textConnection(CVs)
data   <- read.table(txtcon, header = TRUE, sep = "|",
                     strip.white = TRUE, as.is = TRUE)
close(txtcon)
print(CVpooled(data, alpha = 0.20), digits = 4, verbose = TRUE)

## ----PA-----------------------------------------------------------------------
pa.ABE(CV = 0.20, theta0 = 0.92)

## ----DO-----------------------------------------------------------------------
balance <- function(x, y) {
  return(y * (x %/% y + as.logical(x %% y)))
}
do        <- 0.15 # anticipated dropout-rate 15%
seqs      <- 3
n         <- seq(12, 96, 12)
res       <- data.frame(n = n,
                        adj1 = balance(n / (1 - do), seqs), # correct
                        elig1 = NA, diff1 = NA,
                        adj2 = balance(n * (1 + do), seqs), # wrong
                        elig2 = NA, diff2 = NA)
res$elig1 <- floor(res$adj1 * (1 - do))
res$diff1 <- sprintf("%+i", res$elig1 - n)
res$elig2 <- floor(res$adj2 * (1 - do))
res$diff2 <- sprintf("%+i", res$elig2 - n)
invisible(
  ifelse(res$elig1 - n >=0, res$optim <- res$elig1, res$optim <- res$elig2))
res$diff  <- sprintf("%+i", res$optim - n)
names(res)[c(2, 5)] <- c("n'1", "n'2")
res$diff1[which(res$diff1 == "+0")] <- "\u00B10"
res$diff2[which(res$diff2 == "+0")] <- "\u00B10"
res$diff[which(res$diff == "+0")]   <- "\u00B10"
print(res, row.names = FALSE)

## ----lit1---------------------------------------------------------------------
CVfromCI(lower = 0.8323, upper = 1.0392, design = "2x2x4", n = 26)

## ----lit2---------------------------------------------------------------------
n      <- 26
CV.est <- CVfromCI(lower = 0.8323, upper = 1.0392, design = "2x2x4", n = 26)
n.est  <- sampleN.TOST(CV = CV.est, design = "2x2x4", print = FALSE,
                       details = FALSE)[["Sample size"]]
n1     <- balance(seq(n, 18, -1), 2) / 2
n2     <- n - n1
nseqs  <- unique(data.frame(n1 = n1, n2 = n2, n = n))
res    <- data.frame(n1 = nseqs$n1, n2 = nseqs$n2, CV.true = NA,
                     CV.est = CV.est, n.true = NA, n.est = n.est)
for (i in 1:nrow(res)) {
  res$CV.true[i] <- CVfromCI(lower = 0.8323, upper = 1.0392, design = "2x2x4",
                             n = c(res$n1[i], res$n2[i]))
  res$n.true[i]  <- sampleN.TOST(CV = res$CV.true[i], design = "2x2x4",
                                 print = FALSE, details = FALSE)[["Sample size"]]
  res$n.est[i]   <- sampleN.TOST(CV = CV.est, design = "2x2x4",
                                 print = FALSE, details = FALSE)[["Sample size"]]
}
print(signif(res, 5), row.names = FALSE)

## ----dev----------------------------------------------------------------------
CV   <- 0.21
d    <- 0.05 # delta 5%, direction unknown
n    <- sampleN.TOST(CV = CV, theta0 = 1 - d, print = FALSE,
                     details = FALSE)[["Sample size"]]
res1 <- data.frame(CV = CV, theta0 = c(1 - d, 1 / (1 - d)),
                   n = n, power = NA)
for (i in 1:nrow(res1)) {
  res1$power[i] <- power.TOST(CV = CV, theta0 = res1$theta0[i], n = n)
}
n    <- sampleN.TOST(CV = CV, theta0 = 1 + d, print = FALSE,
                     details = FALSE)[["Sample size"]]
res2 <- data.frame(CV = CV, theta0 = c(1 + d, 1 / (1 + d)),
                   n = n, power = NA)
for (i in 1:nrow(res1)) {
  res2$power[i] <- power.TOST(CV = CV, theta0 = res2$theta0[i], n = n)
}
res <- rbind(res1, res2)
print(signif(res[order(res$n, res$theta0), ], 4), row.names = FALSE)

