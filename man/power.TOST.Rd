\name{power.TOST}
\alias{power.TOST}
\title{
Power of the TOST procedure
}
\description{
Calculates the power of the two-one-sided t-tests procedure for various
study designs used in BE studies
}
\usage{
power.TOST(alpha = 0.05, logscale = TRUE, theta1 = 0.8, theta2, diff = 0.95, 
           CV, n, design = "2x2", exact = TRUE)
}
\arguments{
  \item{alpha}{
  Level of significance. Mostly set to 0.05.
  }
  \item{logscale}{
  TRUE or FALSE. Are the evaluation planned after log-transformation?
  }
  \item{theta1}{
  Theta1 is the lower bioequivalence limit.
  In case of logtrans=TRUE it is given as ratio, otherwise as diff. to 1.
  }
  \item{theta2}{
  Theta2 is the upper bioequivalence limit.
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
  \item{n}{
  Number of subjects under study.
  }
  \item{design}{
  Character string describing the study design.
  See known.designs() for designs covered in this package.
  }
  \item{exact}{
  Shall the calculations done exact (based on formulas with Owen's Q) or
  approximately via non-central Student's t-distribution? TRUE or FALSE.
  }
}
\details{
The calculations of power are exact based on Owen's Q-function.
There is also an implementation of approximate power via non-central 
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
"Sample Size Determination for Bioequivalence Assessment \cr 
by Means of Confidence Intervals" \cr 
Int. J. of Clinical Pharmacology, Therapy and Toxicology, 29, 1-8 
}
\author{
D. Labes
}
\seealso{
 \code{\link{sampleN.TOST}, \link{known.designs}}
}
\examples{
# power for the 2x2 cross-over
power.TOST(CV=0.25, n=24)
}
