\name{exppower.TOST}
\alias{exppower.TOST}
\title{
'Expected' power of TOST procedure
}
\description{
Calculates the 'expected' power according to Julious for a variety of study
designs used in bioequivalence studies.
}
\usage{
exppower.TOST(alpha = 0.05, logscale=TRUE, theta0, theta1, theta2,  
              CV, dfCV, n, design = "2x2", robust=FALSE)
}
\arguments{
  \item{alpha}{
Level of significance. Commonly set to 0.05.
  }
  \item{logscale}{
Should the data used on log-transformed or on original scale? TRUE or FALSE.\cr
Defaults to TRUE.
  }  
  \item{theta0}{
'True' or assumed bioequivalence ratio or difference. 
Typically set to 0.95 (default if missing) if \code{logscale=TRUE}.\cr
Defaults to 0.05 if \code{logscale=FALSE}.
  }
  \item{theta1}{
Lower bioequivalence limit as ratio if \code{logscale=TRUE} or as difference.\cr
Can be missing. Defaults then to 0.8 if \code{logscale=TRUE} or to -0.2 if 
\code{logscale=FALSE}.
  }
  \item{theta2}{
Upper bioequivalence limit as ratio if \code{logscale=TRUE} or as difference.\cr
If not given \code{theta2} will be calculated as \code{1/theta1} if \code{logscale=TRUE},
else as \code{-theta1}.
  }
  \item{CV}{
Coefficient of variation as ratio. 
  }
  \item{dfCV}{
Degrees of freedom for the CV (error/residual degree of freedom).
}
  \item{n}{
Number of subjects under study.\cr
Is total number if given as scalar, else number of subjects in the (sequence) 
groups. In the latter case the length of n vector has to be equal to the 
number of (sequence) groups.
}
  \item{design}{
Character string describing the study design.
See \code{known.designs()} for designs covered in this package.
}
  \item{robust}{
Defaults to FALSE.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' evaluation
(aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.
}
}
\details{
This function calculates the so-called 'expected' power based on S.A. Julious
taking into account that usually the CV is not known but estimated from a
previous study / studies with an uncertainty. See references.
}
\value{
Value of expected power according to the input.
}
\references{
S.A. Julious, R.J. Owen\cr
"Sample size calculations for clinical studies allowing for 
uncertainty in variance"\cr
Pharmaceutical Statistics (2006), 5, 29-37

S.A. Julious\cr
"Sample sizes for Clinical Trials"\cr
CRC Press, Chapman & Hall 2010
}
\author{
D. Labes
}
\note{
This function is now also implemented for \code{logscale=FALSE}.\cr
Previous versions only had \code{logscale=TRUE} implemented with lack of this argument.
}

\seealso{
 \code{\link{expsampleN.TOST}, \link{power.TOST}} 
}
\examples{
# expected power for a 2x2 crossover
# CV 30\% known from a pilot study with 12 subjects (-> dfCV=10)
# using all the defaults for other parameters
# should give: [1] 0.735977
exppower.TOST(CV=0.3, dfCV=10, n=40)

# Compare this to the usual power (CV known, "carved in stone")
# gives: [1] 0.8158453
power.TOST(CV=0.3, n=40)
}
