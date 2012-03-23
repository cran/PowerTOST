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
power.TOST(alpha = 0.05, logscale = TRUE, theta1, theta2, theta0, CV, n, 
           design = "2x2", method="exact", robust=FALSE)
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. Conventionally mostly set to 0.05.
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
  \item{CV}{
Coefficient of variation as ratio.\cr
In case of cross-over studies this is the within-subject CV, \cr
in case of a parallel-group design the CV of the total variability.
}
  \item{n}{
Number of subjects under study.\cr 
This is is from package version >0.9-0 the total number of subjects.
}
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
  }
  \item{method}{
Defaults to "exact" in which case the calculation is done based on formulas 
with Owen's Q.\cr
The exact calculation can also be choosen with \code{method="owenq"}
Approximate calculations can be choosen via \code{method="noncentral"} or 
\code{method="nct"} for the approximation using the non-central t-distribution\cr
or via \code{method="central"} or \code{method="shifted"} for the 
approximation via 'shifted' central t-distribution.\cr
The strings for \code{method} may be abbreviated.
  }
  \item{robust}{
Defaults to FALSE. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' evaluation
(aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.\cr
Has only effect for higher-order crossover designs.
  }
}
\details{
The exact calculations of power are based on Owen's Q-function.\cr
The approximate power is implemented via non-central t-distribution\cr
or via 'shifted' central t-distribution.\cr\cr
The formulas used assume balanced studies, i.e. equal number of subjects in
the (sequence) groups.\cr\cr
The formulas for the paired means 'design' do not take a correlation parameter 
into account.\cr 
They are solely based on the paired t-test (TOST of differences = zero).
}
\value{
Value of power according to the input arguments.
}
\section{Warning }{
The function does not vectorize if design is a vector.\cr
The function vectorize properly if n or CV or theta0 are vectors.\cr
Other vector input is not tested yet.
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
# power for the 2x2 cross-over design with 24 subjects
# using all the other default values
# should give: [1] 0.7391155
power.TOST(CV=0.25, n=24)
}
