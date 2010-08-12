\name{power.TOST}
\alias{power.TOST}
\title{
Power of the classical TOST procedure
}
\description{
Calculates the exact or approximate power of the two-one-sided t-tests 
procedure for various study designs used in BE studies.
}
\usage{
power.TOST(alpha = 0.05, logscale = TRUE, theta1, theta2, theta0, diff, CV, 
           n, design = "2x2", exact = TRUE)
}
\arguments{
  \item{alpha}{
Type I error, significance level. Mostly set to 0.05.
}
  \item{logscale}{
Should the data used on log-transformed or on original scale? TRUE or FALSE.\cr
Defaults to TRUE.
}
  \item{theta1}{
Lower bioequivalence limit.\cr
In case of \code{logscale=TRUE} it is given as ratio, otherwise as diff. to 1.\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
}
  \item{theta2}{
Upper bioequivalence limit.\cr
If not given theta2 will be calculated as \code{1/theta1} if \code{logscale=TRUE} \cr
or as \code{-theta1} if \code{logscale=FALSE}.
}
  \item{theta0}{
'True' or assumed bioequivalence ratio. \cr
In case of \code{logscale=TRUE} it must be given as ratio, \cr
otherwise as difference to 1. See examples.\cr
Defaults to 0.95 if \code{logscale=TRUE} or to 0.05 if \code{logscale=FALSE}
}
  \item{diff}{
Same as \code{theta0}. This argument is only for back compatibility with former versions. \cr
It will be ceased in the next version.\cr
If both \code{theta0} and \code{diff} are given, \code{theta0} rules.
}
  \item{CV}{
Coefficient of variation as ratio.\cr
In case of cross-over studies this is the within-subject CV, \cr
in case of a parallel-group design the CV of the total variability.
}
  \item{n}{
Number of subjects under study.
}
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
  }
  \item{exact}{
Shall the calculations done exact (based on formulas with Owen's Q) or
approximately via non-central Student's t-distribution? \cr 
TRUE or FALSE. Defaults to TRUE.
  }
}
\details{
The exact calculations of power are based on Owen's Q-function.\cr
The approximate power is an implementation of  via non-central 
t-distribution.
}
\value{
Value of power according to the input.
}
\references{
Phillips, K. F. (1990) \cr 
"Power of the Two One-Sided Tests Procedure in Bioequivalence" \cr
Journal of Pharmacokinetics and Biopharmaceutics, 18, 137-144.

Diletti D., Hauschke D., and Steinijans V. W. (1991) \cr
"Sample Size Determination for Bioequivalence Assessment by Means of 
Confidence Intervals" \cr 
Int. J. of Clinical Pharmacology, Therapy and Toxicology, 29, 1-8 

See here for a short description:\cr
\url{../doc/BE_power_sample_size_excerpt.pdf}.
}
\author{
D. Labes
}
\seealso{
 \code{\link{sampleN.TOST}, \link{known.designs}}
}
\examples{
# power for the 2x2 cross-over 
# using all the other default values
power.TOST(CV=0.25, n=24)
# gives: [1] 0.7391155
}
