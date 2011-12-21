\name{power2.TOST}
\alias{power2.TOST}
\title{
Power of the classical TOST procedure with unbalanced (sequence) groups
}
\description{
Calculates the exact or approximate power of the two-one-sided t-tests 
procedure for various study designs used in BE studies.
  }
\usage{
power2.TOST(alpha = 0.05, logscale = TRUE, theta1, theta2, theta0, CV, n, 
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
Number of subjects in the (sequence) groups under study.\cr
Must be a vector of length = (sequence) groups.
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
The approximate power is implemented via non-central t-distribution \cr
or via 'shifted' central t-distribution.\cr
This function does'nt make many sense in case of \code{design="paired"}.
}
\value{
Value of power according to the input arguments.
}
\references{
Phillips, K. F. (1990) \cr 
"Power of the Two One-Sided Tests Procedure in Bioequivalence" \cr
Journal of Pharmacokinetics and Biopharmaceutics, 18, 137-144.

Diletti D., Hauschke D., and Steinijans V. W. (1991) \cr
"Sample Size Determination for Bioequivalence Assessment by Means of 
Confidence Intervals" \cr 
Int. J. of Clinical Pharmacology, Therapy and Toxicology, 29, 1-8 

}
\author{
D. Labes
}
\seealso{
 \code{\link{power.TOST}, \link{known.designs}}
}
\examples{
# power for the 2x2 cross-over design with 24 subjects balanced
# using all the other default values
# should give: [1] 0.7391155
power2.TOST(CV=0.25, n=c(12,12))
# the same result 
power.TOST(CV=0.25, n=24)

# power for the 2x2 cross-over design with 24 subjects and
# 2 dropouts in the same sequence group
# should give: [1] 0.6912935
power2.TOST(CV=0.25, n=c(10,12))

}
