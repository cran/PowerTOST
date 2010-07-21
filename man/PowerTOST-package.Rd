\name{PowerTOST-package}
\alias{PowerTOST-package}
\alias{PowerTOST}
\docType{package}
\title{
Power and sample size based on two-one-sided t-test (TOST) procedure
for bioequivalence studies
}
\description{
Contains functions to calculate power and sample size for various study designs 
used for bioequivalence studies.\cr 
See function known.designs() for study designs covered. \cr
Moreover the package contains functions for power and sample size based on 
'expected power' in case of uncertain (estimated) variability.
}
\details{
\tabular{ll}{
Package: \tab PowerTOST\cr
Type: \tab Package\cr
Version: \tab 0.6-2\cr
Date: \tab 2010-07-21\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
Main functions are sampleN.TOST() and power.TOST().\cr
If you prefer sample size based on 'expected' power see the functions
exp.sampleN.TOST() and exp.power.TOST().
}
\author{
D. Labes \cr
Maintainer: D. Labes <DetlewLabes@gmx.de>
}
\references{
Phillips, K. F. (1990)\cr
"Power of the Two One-Sided Tests Procedure in Bioequivalence\cr
Journal of Pharmacokinetics and Biopharmaceutics, 18, 137-144.

Diletti, D., Hauschke, D., and Steinijans, V. W. (1991)\cr
"Sample Size Determination for Bioequivalence Assessment\cr
by Means of Confidence Intervals"\cr
Int. J. of Clinical Pharmacology, Therapy and Toxicology, 29, 1-8

S.A. Joulious\cr
"Sample sizes for Clinical Trials"\cr
CRC Press, Chapman & Hall 2010
}
\keyword{ package }
