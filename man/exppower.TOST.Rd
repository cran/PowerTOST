\name{exppower.TOST}
\alias{exppower.TOST}
\title{
'Expected' power of TOST procedure
}
\description{
Calculates the expected power according to Joulious for a variety of study
designs used in bioequivalence studies
}
\usage{
exppower.TOST(alpha = 0.05, theta1 = 0.8, theta2, diff = 0.95, 
               CV, dfCV, n, design = "2x2")
}
\arguments{
  \item{alpha}{
Level of significance. Mostly set to 0.05.
}
  \item{theta1}{
Theta1 is the lower bioequivalence limit as ratio.
}
  \item{theta2}{
Theta2 is the upper bioequivalence limit.
If not given theta2 will be calculated as 1/theta1.
}
  \item{diff}{
'True' or assumed bioequivalence ratio.
}
  \item{CV}{
Coefficient of variation as ratio. 
}
  \item{dfCV}{
Degrees of freedom for the CV (error/residual degree of freedom).
}
  \item{n}{
Number of subjects to be planned.
}
  \item{design}{
Character string describing the study design.
See known.designs() for designs covered in this package.
}
}
\details{
This function calculates the so-called 'expected' power based on S.A. Joulious
taking into account that usually the CV is not known but estimated from
previous study / studies with an uncertainty. See reference.
}
\value{
Value of power according to the input.
}
\references{
S.A. Joulious, R.J. Owen\cr
"Sample size calculations for clinical studies allowing for \cr
uncertainty in variance"\cr
Pharmaceutical Statistics (2006), 5, 29-37

S.A. Joulious\cr
"Sample sizes for Clinical Trials"\cr
CRC Press, Chapman & Hall 2010
}
\author{
D. Labes
}
\note{
This function is only implemented for logscale=TRUE. Therefore this argument is
not included in the function call.
}

\seealso{
 \code{\link{power.TOST}, \link{expsampleN.TOST}}
}
\examples{
# expected power for a 2x2 crossover
# CV 30\% known from a pilot study with 12 subjects (dfCV=10)
# should give: [1] 0.3575385
exppower.TOST(CV=0.3,dfCV=10,n=20)
}
