# --------------------------------------------------------------------------
# shadowtext() lifted from the package TeachingDemos
# author Greg Snow, orphaned on CRAN 2024-02-10
# request from Prof Brian Ripley to comply with the CRAN policy
# Later on the orphan status was abandoned, but decided to use  
# the recoded function nevertheless
shtxt <- function(x, y = NULL, labels, col = "white", bg = "black",
                  theta = seq(pi / 32, 2 * pi, length.out = 64), r = 0.1, 
                  cex = 1, ... ) 
{
  xy  <- xy.coords(x, y)
  fx  <- grconvertX(xy$x, to = "nfc")
  fy  <- grconvertY(xy$y, to = "nfc")
  fx0 <- r * strwidth("A", units = "figure", cex = cex)
  fy0 <- r * strheight("A", units = "figure", cex = cex)
  for (step in theta) {
    text(grconvertX(fx + cos(step) * fx0, from = "nfc"),
         grconvertY(fy + sin(step) * fy0, from = "nfc"),
         labels, cex = cex, col = bg, ...)
  }
  text(xy$x, xy$y, labels, cex = cex, col = col, ... )
}

# S3 method for printing the results of pa.ABE(), pa.scABE()
# --------------------------------------------------------------------------
print.pwrA <- function(x, digits=4, plotit=TRUE, ...)
{

  if (interactive() && plotit) plot(x)

  min.pwr  <- x$minpower
  CV.max   <- max(x$paCV[, "CV"])
  CVmaxI   <- which(x$paCV[, "CV"]==CV.max)
  CV.min   <- min(x$paCV[, "CV"])
  CVminI   <- which(x$paCV[, "CV"]==CV.min)
  min.pwrN <- min(x$paN[, "pwr"])
  min.N    <- min(x$paN[, "N"])
  incr     <- x$incr

  if (abs(x$paCV[CVmaxI, "pwr"]-min.pwr)>1e-4) CV.max <- NA
  if (abs(x$paCV[CVminI, "pwr"]-min.pwr)>1e-4) CV.min <- NA

  if (x$plan[1, "theta0"]<=1) {
    min.theta0 <- min(x$paGMR[, "theta0"])
  } else {
    min.theta0 <- max(x$paGMR[, "theta0"])
  }
  method <- x$method
  if (method == "scABE") {
    meth <- switch(x$regulator,
                   EMA=   " (EMA/ABEL)",
                   HC =   " (HC/ABEL)",
                   FDA=   " (FDA/RSABE)",
                   GCC=   " (GCC/ABEL)"
                  )
    method <- paste0(method, meth)
  }
  cat("Sample size plan ", method, "\n",sep="")
  # without "nlast"
  print(x$plan[,-ncol(x$plan)], row.names = FALSE)
  cat("\nPower analysis\n")
  cat("CV, theta0 and number of subjects leading to min. acceptable ",
      "power of ~", round(min.pwr, digits), ":\n", sep="")
  #react to NTID where there may be a CV.min, CV.max which
  if (method!="NTID") {
    # let to power=minpower
    cat(" CV= ", round(CV.max, digits), ", theta0= ",
        round(min.theta0, digits),"\n", sep="")
  } else {
      # let to power=minpower
      cat(" CV = (", round(CV.min, digits), ", ", round(CV.max, digits),
          "), theta0= ", round(min.theta0, digits),"\n", sep="")
  }
  cat(" n = ", min.N, " (power= ", round(min.pwrN, digits), ")\n", sep="")
  if (incr) {
    cat("Sample size was increased from the estimated one in order to comply",
        "\nwith regulatory requirements.")
  }
  cat("\n")
}

# ----------------------------------------------------------------------------
# S3 method for plotting the results of pa.ABE(), pa.scABE(), pa.NTIDFDA()
#
# Author D.Labes, from original code by H. Schuetz
# reworked by H. Schuetz to avoid overlay of target power line with legend
# ----------------------------------------------------------------------------
plot.pwrA <- function(x, pct=TRUE, ratiolabel="theta0", cols=c("blue", "red"), ...)
{
  # make colors between the both given
  mkColors <- function(){
    # all varibles are seen if this function is inside plot
    colno <- cut(pwr, breaks=unique(c(1*fact, targetpower, pwr[pwr<=targetpower])),
                 labels=FALSE, include.lowest=TRUE, right=FALSE)
    clr1 <-colorRampPalette(rev(cols))(length(unique(colno)))
    clr  <- clr1[colno]
    return(clr)
  }

  if (pct) {
    fact    <- 100
    ylabtxt <- "power (%)"
    dec     <- 2
    pctsign <- "%"
    GMRmain <- sprintf("%.4g%%", 100*x$plan[1,"theta0"])
  }   else {
    fact    <- 1
    ylabtxt <- "power"
    dec     <- 4
    pctsign <- ""
    GMRmain <- sprintf("%.4g", x$plan[1,"theta0"])
  }

  n <- nrow(x$paCV)
  # Attention next functiones only if last value is minpower
  minpower    <- fact*x$minpower
  CV.max      <- fact*max(x$paCV[, "CV"])
  CV.min      <- fact*min(x$paCV[, "CV"])

  # CV of call of pa.XXX() function
  if (x$method=="ABE") {
    CV  <- fact*x$plan[1,"CV"]
  } else {
    # Attention this functions only if CVwT==CVwR!!!
    CV  <- fact*x$plan[1,"CVwR"]
  }
  GMR         <- x$plan[1, "theta0"]
  theta1      <- x$plan[1, "theta1"]
  theta2      <- x$plan[1, "theta2"]
  n.est       <- x$plan[1, "Sample size"]
  pwr.est     <- fact*x$plan[1, "Achieved power"]
  targetpower <- fact*x$plan[1, "Target power"]
  reg         <- x$regulator
  design      <- x$plan[1, "Design"]

  mklegend <- function(method){
    if (method=="scABE"){
      algo <- ifelse(reg == "FDA", "RSABE", "ABEL")
      legend("topright", legend=c(algo, paste0("(",reg,")")), cex=0.90,
      bg="white", box.lty=0, x.intersp=0)
    }
    if (method=="NTID"){
      legend("topright", legend=c("RSABE", "(NTID)"), cex=0.90,
      bg="white", box.lty=0, x.intersp=0)
    }
  }

  op <- par(no.readonly=TRUE) # save par() options
  on.exit(par(op))
  par(mar =c(c(3.5, 3.5, 2.5, 0.75))+0.1) # bottom, left, top, right
  par(cex.main=0.95, font.main=1, cex.axis=0.95, cex.lab=0.95,
      mgp =c(2, 0.75, 0), tcl = -0.2)

  # plot on a panel of 4 pieces
  split.screen(c(2, 2))

  screen(1) ### 'Sensitivity' of CV (GMR and n constant) ###
  pwr <- as.numeric(fact*x$paCV[,"pwr"])
  CVs <- as.numeric(fact*x$paCV[,"CV"])
  seg <- length(pwr); s <- seq(seg-1)
  clr <- mkColors()
  xlabtxt <- "CV"
  if (fact==100) xlabtxt <- "CV (%)"

  if (x$method == "ABE") {
    plot(CVs, pwr, type="n",
         main=paste0("Higher variability\n",
                     "constant: ", ratiolabel, " = ", GMRmain, ", n = ", n.est),
         lwd=2, xlab=xlabtxt, ylab="", las=1)
    mtext(side=2, ylabtxt, line=2.5)
    grid()
    abline(h=c(targetpower, fact*0.8, minpower), lty=3)
    segments(CVs[s], pwr[s], CVs[s+1], pwr[s+1], lwd=2, col=clr[s])
    points(CVs[1], pwr[1], col=clr[1], pch=16, cex=1.25)
    points(CVs[seg], pwr[seg], col=clr[seg], pch=16, cex=1.25)
    shtxt(CV, (minpower+(pwr.est-minpower)*0.1),
               labels=paste0("CV = ", signif(CV.max, 4), pctsign," (",
                             round(minpower, dec), pctsign,")"),
               col="black", bg="white", pos=4, r=0.5, cex=0.9)

  } else {
    # any scABE (including RSABE NTID)
    plot(CVs, pwr, type="n",
         main=paste0("Lower/higher variability\n",
                     "constant: ", ratiolabel, " = ", GMRmain, ", n = ", n.est),
         lwd=2, xlab=xlabtxt, ylab="", las=1)
    grid()
    abline(h=c(targetpower, 0.8*fact, minpower), lty=3)
    mtext(side=2, ylabtxt, line=2.5)
    mklegend(x$method)
    segments(CVs[s], pwr[s], CVs[s+1], pwr[s+1], lwd=2, col=clr[s])
    # mark the plan CV and power
    points(CV, pwr.est, col=cols[1], pch=16, cex=1.25)
    # mark the max. CV if pwr for it is = minpower
    if (abs(pwr[seg]-minpower)/minpower<=1e-4){
      points(CVs[seg], pwr[seg], col=clr[seg], pch=16, cex=1.25)
    }
    txt <- paste0("CV = ", signif(CV.max, 4), pctsign, " (",
                  round(minpower, dec), pctsign, ")")
    if  (x$method=="NTID") {
      if (abs(pwr[1]-minpower)/minpower <= 1e-4) {
        # we have also CV.min with power=minpower
        points(CV.min, pwr[1], col=clr[seg], pch=16, cex=1.1)
        txt <- paste0("CV = (", signif(CV.min, 4),", ", signif(CV.max, 4),
                      pctsign, ") (",round(minpower, dec), pctsign, ")")
      }
    }
    shtxt(min(CVs), (minpower+(max(pwr)-minpower)*0.1),
               labels=txt,col="black", bg="white", pos=4, r=0.5, cex=0.9)
  }
  box()

  screen(2) ### 'Sensitivity' of GMR (CV and n constant) ###
  pwr <- as.numeric(fact*x$paGMR[, "pwr"])
  GMRs <- as.numeric(x$paGMR[, "theta0"])
  #GMR.min <- ifelse(GMR<=1, GMRs[1], GMRs[length(GMRs)])
  GMR.min <- GMRs[1] # or better min(GMRs) or max(GMRs)?
  seg <- length(pwr); s <- seq(seg-1)
  clr  <- mkColors()
  if (fact == 1) { # like in previous versions (ratio)
    plot(GMRs, pwr, type="n",
         main=paste0("Larger deviation from 1\n",
                     "constant: CV = ", CV, pctsign,", n = ", n.est),
         lwd=2, xlim=c(GMR, GMR.min), xlab=ratiolabel, ylab="", las=1,
         cex.main=0.95, cex.axis=0.95)
    grid()
    abline(h=c(targetpower, fact*0.8, minpower), lty=3)
    mklegend(x$method)
    mtext(ylabtxt, side=2, line=2.5)
    segments(GMRs[s], pwr[s], GMRs[s+1], pwr[s+1], lwd=2, col=clr[s])
    # the next assumes that the values start at GMR and end on GMR.min (maybe also max!)
    # TODO rework if plan.GMR not at border
    points(GMRs[1], pwr[1], col=clr[1], pch=16, cex=1.25)
    points(GMRs[seg], pwr[seg], col=clr[seg], pch=16, cex=1.25)
    shtxt(GMR, (minpower+(pwr.est-minpower)*0.1),
               labels=paste0(ratiolabel, " = ", signif(GMR.min, 4), " (",
                             round(minpower, dec), pctsign, ")"),
               col="black", bg="white", pos=4, r=0.5, cex=0.9)
    } else {
    plot(100*GMRs, pwr, type="n",
         main=paste0("Larger deviation from 100%\n",
                     "constant: CV = ", CV, pctsign,", n = ", n.est),
         lwd=2, xlim=100*c(GMR, GMR.min), xlab=paste(ratiolabel, "(%)"), ylab="", las=1,
         cex.main=0.95, cex.axis=0.95)
    grid()
    abline(h=c(targetpower, fact*0.8, minpower), lty=3)
    mklegend(x$method)
    mtext(ylabtxt, side=2, line=2.5)
    segments(100*GMRs[s], pwr[s], 100*GMRs[s+1], pwr[s+1], lwd=2, col=clr[s])
    # the next assumes that the values start at GMR and end on GMR.min (maybe also max!)
    # TODO rework if plan.GMR not at border
    points(100*GMRs[1], pwr[1], col=clr[1], pch=16, cex=1.25)
    points(100*GMRs[seg], pwr[seg], col=clr[seg], pch=16, cex=1.25)
    text(100*GMR, (minpower+(pwr.est-minpower)*0.1),
         labels=paste0(ratiolabel, " = ",signif(100*GMR.min, 4), "% (",
                       round(minpower, dec), pctsign, ")"),
         cex=0.9, pos=4)
    shtxt(100*GMR, (minpower+(pwr.est-minpower)*0.1),
               labels=paste0(ratiolabel, " = ", signif(100*GMR.min, 4), "% (",
                             round(minpower, dec), pctsign, ")"),
               col="black", bg="white", pos=4, r=0.5, cex=0.9)
  }
  box()

  screen(3) ### Sensitivity of n (GMR and CV constant) ###
  pwr <- as.numeric(fact*x$paN[, "pwr"])
  Ns  <- as.numeric(x$paN[, "N"])
  clr <- mkColors()
  if (length(clr) == 1) clr <- cols[1]
  xticks <- NULL
  nNs    <- length(Ns)
  if (nNs<5 & nNs>1) xticks <- c(max(Ns), min(Ns), nNs-1)
  plot(Ns, pwr, type="n",
       main=paste0("Dropouts\n",
                   "constant: ", ratiolabel, " = ", GMRmain,
                   ", CV = ", CV, pctsign),
       lwd=2, xlim=c(max(Ns), min(Ns)), ylim=c(minpower, pwr.est),
       xlab="n", xaxp=xticks,
       ylab="", las=1, cex.main=0.95)
  grid()
  abline(h=c(targetpower, fact*0.8, minpower), lty=3)
  mklegend(x$method)
  mtext(side=2, ylabtxt, line=2.5)
  points(Ns, pwr, pch=16, cex=0.8, col=clr)
  points(Ns[length(Ns)], pwr[length(Ns)], col=clr[length(Ns)],
         pch=16, cex=1.25)
  points(n.est, pwr.est, col=clr[1], pch=16, cex=1.25)
  # label even drop-outs
  do.even  <- Ns[c(TRUE, FALSE)]
  pwr.even <- pwr[c(TRUE, FALSE)]
  shtxt(do.even, pwr.even, labels=n.est-do.even,
             col="black", bg="white", pos=3, r=0.4, cex=0.75)
  shtxt(max(Ns), (minpower+(pwr.est-minpower)*0.1),
             labels=paste0("n = ", min(Ns), " (", signif(min(pwr), 4), pctsign, ")"),
             col="black", bg="white", pos=4, r=0.5, cex=0.9)
  box()

  screen(4) ### Some basic information ###
  if (x$method != "NTID") {
    if (fact == 1) {
      CVtxt <- sprintf("  CV = %.4f (%+5.1f%%)",
                     CV.max, 100*(CV.max-CV)/CV)
    } else {
      CVtxt <- sprintf("  CV = %5.2f%% (%+5.1f%%)",
                     CV.max, 100*(CV.max-CV)/CV)
    }
  } else {
    # CVtxt <- "" # why?
    if (abs(fact*x$paCV[1, "pwr"]-minpower)/minpower <= 1e-4) {
      # we have also CV.min with power=minpower
      if (fact == 1) {
        CVtxt <- sprintf("  CVmin = %.5f (%+5.1f%%)", CV.min, 100*(CV.min-CV)/CV)
        CVtxt <- c(CVtxt, sprintf("  CVmax = %.4f (%+5.1f%%)", CV.max, 100*(CV.max-CV)/CV))
      } else {
        CVtxt <- sprintf("  CVmin = %5.3f%% (%+5.1f%%)", CV.min, 100*(CV.min-CV)/CV)
        CVtxt <- c(CVtxt, sprintf("  CVmax = %5.2f%% (%+5.1f%%)", CV.max, 100*(CV.max-CV)/CV))
      }
    } else {
      # we have only CV.max with power=minpower
      if (fact == 1) {
        CVtxt <- sprintf("  CV = %.4f (%+5.1f%%)", CV.max, 100*(CV.max-CV)/CV)
      } else {
        CVtxt <- sprintf("  CV = %5.2f%% (%+5.1f%%)", CV.max, 100*(CV.max-CV)/CV)
      }
    }
  }
 if (x$method=="ABE") {
   BEARtxt <- "  BE margins:"
   if (fact == 1) { # ratios
     BEARtxt <- c(BEARtxt, sprintf("    %.4f %s %.4f",
                                  theta1, "...", theta2))
   } else { # percent
     BEARtxt <- c(BEARtxt, sprintf("    %.2f%% %s %.2f%%",
                                  100*theta1, "...", 100*theta2))
   }
 }
 if (x$method == "scABE") {
   # (widened) acceptance range
   if (x$regulator == "FDA"){
     Ltxt <-"  implied BE margins:"
     wtheta1 <- min(theta1,exp(CV2se(CV/fact)*log(theta1)/0.25))
     wtheta2 <- max(theta2,exp(CV2se(CV/fact)*log(theta2)/0.25))
   } else { #EMA
     Ltxt <- "  (widened) BE margins:"
     CVV <- min(0.5,CV/fact)      # cap
     wtheta1 <- min(theta1,exp(-CV2se(CVV)*0.76))
     wtheta2 <- max(theta2,exp(CV2se(CVV)*0.76))
   }
   if (fact == 1) { # ratios
     BEARtxt <- c(Ltxt, sprintf("    %.4f %s %.4f",
                                wtheta1, "...", wtheta2))
   } else { # percent
     BEARtxt <- c(Ltxt, sprintf("    %.2f%% %s %.2f%%",
                                100*wtheta1, "...", 100*wtheta2))
   }
 }
 if (x$method == "NTID") {
   Ltxt <-"  implied BE margins:"
   wtheta1 <- max(theta1,exp(CV2se(CV/fact)*log(0.9)/0.1))
   wtheta2 <- min(theta2,exp(-CV2se(CV/fact)*log(0.9)/0.1))
   if (fact == 1) { # ratios
     BEARtxt <- c(Ltxt, sprintf("    %.4f %s %.4f",
                                wtheta1, "...", wtheta2))
   } else { # percent
     BEARtxt <- c(Ltxt, sprintf("    %.2f%% %s %.2f%%",
                                100*wtheta1, "...", 100*wtheta2))
   }
 }
 plot(1, type="n", axes=F, xlab="", ylab="")
 if (x$method != "NTID") {
   cex <- 0.9
 } else {
   cex <- 0.85 # more lines require smaller font
 }
 if (fact == 100) { # percent
   legend("topleft", inset=-0.055,
          legend=c(paste0(design, " design", "; assumed:"),
                   sprintf("  %s %.2f%%%s%s%s %.2f%%", "CV =", CV, ", ", ratiolabel, " =", 100*GMR),
                   BEARtxt,
                   "power:",
                   sprintf("  %s %2.0f%%", "target =", targetpower),
                   sprintf("  %s %5.2f%% %s %i%s", "estimated =", pwr.est,
                           "(n =", n.est, ")"),
                   sprintf("  %s %2.0f%%", "min. acceptable =", minpower),
                   "acceptable (relative) deviations:",
                   #TODO:react to RSABE NTID where there may be also a CVmin
                   CVtxt,
                   sprintf("  %s%s %5.2f%% (%+5.2f%%)",
                           ratiolabel, " =", 100*GMR.min, 100*(GMR.min-GMR)/GMR),
                   sprintf("  %s %i (%+5.1f%%)",
                           "n =", min(Ns), 100*(min(Ns)-n.est)/n.est)),
          bty="n", cex=cex)
  } else { # ratios
    legend("topleft", inset=-0.055,
           legend=c(paste0(design, " design", "; assumed:"),
                    sprintf("  %s %.4f%s%s%s %.4f", "CV =", CV, ", ", ratiolabel, " =", GMR),
                    BEARtxt,
                    "power:",
                    sprintf("  %s %5.4f", "target =", targetpower),
                    sprintf("  %s %5.4f %s %i%s", "estimated =", pwr.est,
                            "(n =", n.est, ")"),
                    sprintf("  %s %5.4f", "min. acceptable =", minpower),
                    "acceptable (relative) deviations:",
                    CVtxt,
                    sprintf("  %s%s %.4f (%+5.2f%%)",
                            ratiolabel, " =", GMR.min, 100*(GMR.min-GMR)/GMR),
                    sprintf("  %s %i (%+5.1f%%)",
                            "n =", min(Ns), 100*(min(Ns)-n.est)/n.est)),
           bty="n", cex=cex)
  }

  close.screen(all.screens=TRUE)
}
