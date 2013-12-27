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
(based on Fieller's confidence (fiducial) interval).\cr 
These functions are intended for studies with clinical endpoints.\cr\cr
Contains further functions for power and sample size calculations based on 
non-inferiority test. \cr 
This is not a TOST procedure (but rather OOST ;-)) but eventually useful if the 
question of 'non-superiority' must be evaluated within a BE study.\cr
The power and sample size calculations based on non-inferiority test may also 
performed via 'expected' power in case of uncertain (estimated) variability.\cr\cr
Contains functions to calculate power and sample size for the BE decision 
via scaled (widened) BE acceptance limits based on simulations.\cr
Contains functions to calculate power and sample size for the 
BE decision via linearized reference scaled ABE criterion according to the FDA 
favored procedure based on simulations.\cr
Contains further functions to calculate power and sample size according to the 
FDA procedure for NTID's, also via simulations.
}
\details{
\tabular{ll}{
Package: \tab PowerTOST\cr
Type: \tab Package\cr
Version: \tab 1.1-08\cr
Date: \tab 2013-12-27\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes
}
Main functions are \code{sampleN.TOST()} and \code{power.TOST()} for usual 
power and sample size calculations.\cr
If you prefer sample size based on 'expected' power see the functions
\code{expsampleN.TOST()} and \code{exppower.TOST()}.\cr\cr
The main functions for equivalence of the ratio of means with normality on the
original scale are \code{power.RatioF()} and \code{sampleN.RatioF()}.\cr\cr
The functions for calculating power and sample size for the non-inferiority case
are \code{power.noninf()} and \code{sampleN.noninf()}.\cr\cr
The functions for calculating 'expected' power and sample size for the 
non-inferiority case are \code{exppower.noninf()} and \code{expsampleN.noninf()}.\cr\cr
The functions for calculating 'empirical' power and sample size for scaled
average bioequivalence are \code{power.scABEL()} and \code{sampleN.scABEL()}
for the EMA adopted method via widened BE acceptance limits or 
\code{power.RSABE()} and \code{sampleN.RSABE()} for the FDA favored method via
linearized scaled ABE criterion, respectively. These functions are based on
simulations.\cr
The functions for calculating 'empirical' power and sample size via the
FDA method for NTID's are \code{power.NTIDFDA()} and \code{sampleN.NTIDFDA()}. 
These functions are based on simulations.\cr\cr
The package contains further some utility functions (see Index).
}
\author{
D. Labes \cr
with inspirations from Kem Phillips function \code{power.equivalence.md()} 
in package MBESS.

Maintainer: Detlew Labes <DetlewLabes at gmx.de>
}
\references{
Phillips, K. F. (1990)\cr
"Power of the Two One-Sided Tests Procedure in Bioequivalence"\cr
Journal of Pharmacokinetics and Biopharmaceutics, 18, 137-144.

Diletti, D., Hauschke, D., and Steinijans, V. W. (1991)\cr
"Sample Size Determination for Bioequivalence Assessment\cr
by Means of Confidence Intervals"\cr
Int. J. of Clinical Pharmacology, Therapy and Toxicology, 29, 1-8

S.A. Julious, R.J. Owen (2006)\cr
"Sample size calculations for clinical studies allowing for \cr
uncertainty in variance"\cr
Pharmaceutical Statistics (2006), 5, 29-37

S.A. Julious (2010)\cr
"Sample sizes for Clinical Trials"\cr
CRC Press, Chapman & Hall 2010

Hauschke D., Kieser M., Diletti E. and Burke M. (1999)\cr
"Sample size determination for proving equivalence based on the ratio \cr
of two means for normally distributed data" \cr
Stat. Med. 18(1) p93-105 (1999) 

Hauschke D., Steinijans V. and Pigeot I.\cr
"Bioequivalence studies in Drug Development"\cr
John Wiley & Sons, Chichester (2007)\cr
Chapter 5 and 10.3

FDA "Draft Guidance on Progesterone"\cr
Recommended Apr 2010; Revised Feb 2011\cr
http://www.fda.gov/downloads/Drugs/GuidanceComplianceRegulatoryInformation/Guidances/UCM209294.pdf

Laszlo Tothfalusi and Laszlo Endrenyi\cr
"Sample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs"\cr
J. Pharm. Pharmaceut. Sci. (www.cspsCanada.org) 15(1) 73 - 84, 2011

Tothfalusi L., Endrenyi L. and A. Garcia Arieta\cr
"Evaluation of Bioequivalence for Highly Variable Drugs with Scaled Average 
Bioequivalence"\cr
Clin. Pharmacokin. 48/11, 725-743 (2009)

FDA "Draft Guidance on Warfarin Sodium"\cr
Recommended Dec 2012\cr
http://www.fda.gov/downloads/Drugs/GuidanceComplianceRegulatoryInformation/Guidances/UCM201283.pdf

BEBAC forum: \cr
categories 'Power/Sample size' and 'R for BE/BA'\cr
\url{http://forum.bebac.at}
}
\keyword{ package }
