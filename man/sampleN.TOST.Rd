\name{sampleN.TOST}
\alias{sampleN.TOST}
\title{
Sample size based on power of TOST
}
\description{
Calculates the necessary sample size to have at least a given power.
}
\usage{
sampleN.TOST(alpha = 0.05, targetpower = 0.8, logscale = TRUE, 
             theta1 = 0.8, theta2, diff = 0.95, CV, design = "2x2", 
             exact = TRUE, print = TRUE, details = FALSE)
}
\arguments{
  \item{alpha}{
Error probability. Typically set to 0.05.
}
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1.
Typical values are 0.8 or 0.9.
}
  \item{logscale}{
Should the data used on log-transformed or on original scale? TRUE or FALSE.
Defaults to TRUE.
}
  \item{theta1}{
Lower bioequivalence limit.
In case of logscale=TRUE it must be given as ratio, 
otherwise as difference to 1.
}
  \item{theta2}{
Upper bioequivalence limit.
If not given theta2 will be calculated as 1/theta1 if logscale=TRUE or
as -theta1 if logscale=FALSE.
}
  \item{diff}{
'True' or assumed bioequivalence ratio. In case of logscale=TRUE it must be 
given as ratio, otherwise as difference to 1. See examples.
}
  \item{CV}{
Variation coefficient as ratio.
}
  \item{design}{
Character string describing the study design.\cr
See known.designs() for designs covered in this package.
}
  \item{exact}{
Shall the calculations done exact (based on formulas with Owen's Q) \cr
or approximately via non-central Student's t-distribution? TRUE or FALSE. \cr
Defaults to TRUE.
}
  \item{print}{
If TRUE (default) the function prints its results. If FALSE only
a data.frame with the results will be returned. 
}
  \item{details}{
If TRUE the design characteristics and the steps during
sample size calculations will be shown. \cr Defaults to FALSE.
}
}
\value{
A data.frame with the results will be returned.
}
\references{
Phillips, K. F. (1990)\cr 
"Power of the Two One-Sided Tests Procedure in Bioequivalence"\cr
Journal of Pharmacokinetics and Biopharmaceutics, 18, 137-144.

Diletti, D., Hauschke, D., and Steinijans, V. W. (1991)\cr
"Sample Size Determination for Bioequivalence Assessment\cr
by Means of Confidence Intervals"\cr
Int. J. of Clinical Pharmacology, Therapy and Toxicology, 29, 1-8. 
}
\author{
D. Labes
}
\seealso{
  \code{\link{power.TOST}, \link{known.designs}}
}
\examples{
# Exact calculation for a classical 2x2 cross-over
# BE limits 80 ... 125\%, assumed true BE ratio 0.95
# intra-subject CV=30\%
# using all the default values
sampleN.TOST(CV=0.3)

# Exact calculation for a parallel group design
# evaluation on the original (untransformed) scale
# BE limits 80 ... 120\%, assumed true BE ratio 0.95
# total CV=20\%
sampleN.TOST(CV=0.2, logscale=FALSE, theta1=-0.2, diff=-0.05, design="parallel")
}
