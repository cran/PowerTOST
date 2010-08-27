\name{PowerTOST-package}
\alias{PowerTOST-package}
\alias{PowerTOST}
\docType{package}
\title{
Power and sample size based on two one-sided t-test (TOST) procedure
for bioequivalence studies
}
\description{
Contains functions to calculate power and sample size for various study designs 
used for bioequivalence studies.\cr 
See function \code{known.designs()} for study designs covered. \cr\cr
Moreover the package contains functions for power and sample size based on 
'expected power' in case of uncertain (estimated) variability.\cr\cr
Added are functions for the power and sample size for the ratio of 
two means with normally distributed data on the original scale 
(based on Fieller's confidence (fiducial) interval).\cr These functions are
intended for studies with clinical endpoints.
}
\details{
\tabular{ll}{
Package: \tab PowerTOST\cr
Type: \tab Package\cr
Version: \tab 0.7-2\cr
Date: \tab 2010-08-27\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
Main functions are \code{sampleN.TOST()} and \code{power.TOST()} for usual 
power and sample size calculations.\cr
If you prefer sample size based on 'expected' power see the functions
\code{expsampleN.TOST()} and \code{exppower.TOST()}.\cr
The main functions for equivalence of the ratio of means with normality on the
original scale are \code{power.RatioF()} and \code{sampleN.RatioF()}.
}
\author{
D. Labes \cr
with inspirations from Kem Phillips function \code{power.equivalence.md()} 
in package MBESS.

Maintainer: D. Labes <DetlewLabes@gmx.de>
}
\references{
Phillips, K. F. (1990)\cr
"Power of the Two One-Sided Tests Procedure in Bioequivalence"\cr
Journal of Pharmacokinetics and Biopharmaceutics, 18, 137-144.

Diletti, D., Hauschke, D., and Steinijans, V. W. (1991)\cr
"Sample Size Determination for Bioequivalence Assessment\cr
by Means of Confidence Intervals"\cr
Int. J. of Clinical Pharmacology, Therapy and Toxicology, 29, 1-8

S.A. Joulious, R.J. Owen (2006)\cr
"Sample size calculations for clinical studies allowing for \cr
uncertainty in variance"\cr
Pharmaceutical Statistics (2006), 5, 29-37

S.A. Joulious (2010)\cr
"Sample sizes for Clinical Trials"\cr
CRC Press, Chapman & Hall 2010

Hauschke D., Kieser M., Diletti E. and Burke M. (1999)\cr
"Sample size determination for proving equivalence based on the ratio \cr
of two means for normally distributed data" \cr
Stat. Med. 18(1) p93-105 (1999) 

Hauschke D., Steinijans V. and Pigeot I.\cr
"Bioequivalence studies in Drug Development"\cr
Chapter 5 and 10.3, John Wiley & Sons, Chichester (2007)
}
\keyword{ package }
