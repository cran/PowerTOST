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
exppower.TOST(alpha = 0.05, theta1 = 0.8, theta2, theta0 = 0.95, 
              CV, dfCV, n, design = "2x2")
}
\arguments{
  \item{alpha}{
Level of significance. Commonly set to 0.05.
}
  \item{theta1}{
Lower bioequivalence limit as ratio.
}
  \item{theta2}{
Upper bioequivalence limit as ratio.
If not given \code{theta2} will be calculated as \code{1/theta1}.
}
  \item{theta0}{
'True' or assumed bioequivalence ratio. Tpically set to 0.95.
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
See \code{known.designs()} for designs covered in this package.
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
This function is only implemented for logscale=TRUE.\cr
Therefore this argument is not included in the function call.
}

\seealso{
 \code{\link{power.TOST}, \link{expsampleN.TOST}}
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
