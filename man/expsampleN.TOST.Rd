\name{expsampleN.TOST}
\alias{expsampleN.TOST}
\title{
Sample size based on 'expected' power
}
\description{
Calculates the sample size based on Joulious 'expected' power for a
variety of study designs used in bioequivalence studies. See known.designs() 
for the study designs covered
}
\usage{
expsampleN.TOST(alpha = 0.05, targetpower = 0.8, theta1 = 0.8, theta2, 
                theta0 = 0.95, CV, dfCV, alpha2 = 0.05, design = "2x2", 
                print = TRUE, details = FALSE)
}
\arguments{
  \item{alpha}{
Error probability. Typically set to 0.05.
}
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1.
Typical values are 0.8 or 0.9.
}
  \item{theta1}{
Lower bioequivalence limit as ratio.
}
  \item{theta2}{
Upper bioequivalence limit as ratio.
If not given \code{theta2} will be calculated as \code{1/theta1}.
}
  \item{theta0}{
'True' or assumed bioequivalence ratio. Typically set to 0.95.
}
  \item{CV}{
Coefficient of variation as ratio. \cr
May be given as vector. Then the CV's were pooled as weighted mean with
their df=degrees of freedom as weights.
}
  \item{dfCV}{
Degrees of freedom for the CV's. Must be a vector of same length as CV.\cr
The df's are usually (in case of cross-over studies) \code{N-nseq} where N is the 
total sample size and nseq the number of sequences.
}
  \item{alpha2}{
Significance level for the calculation of an upper confidence interval
for CV (for informal purposes only).
}
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
}
  \item{print}{
If TRUE (default) the function prints its results.\cr
If FALSE only a data.frame with the results will be returned.
}
  \item{details}{
If TRUE the design characteristics and the steps during
sample size calculations will be shown.\cr Defaults to FALSE.
}
}
\details{
The sample size is calculated based on 'expected' power via Joulious formulas 
based on non-central t-distribution.
}
\value{
A data.frame with the input values and the result of the sample size.
}
\references{
S.A. Joulious, R.J. Owen\cr
"Sample size calculations for clinical studies allowing for 
uncertainty in variance"\cr
Pharmaceutical Statistics (2006), 5, 29-37

S.A. Joulious\cr
"Sample sizes for Clinical Trials"\cr
CRC Press, Chapman & Hall, Boca Raton 2010
}
\author{
D. Labes
}
\note{
The returned sample size is in case of a parallel-group design the number of
subjects in each treatment group.\cr
This function is only implemented for logscale=TRUE. Therefore this argument is
not included in the function definition.
}

\seealso{
 \code{\link{exppower.TOST}, \link{known.designs}, \link{sampleN.TOST}}
}
\examples{
# Classical 2x2 cross-over, target power = 80\%,
# BE limits 80 ... 125\%, assumed true BE ratio = 95%,
# intra-subject CV=30\% estimated with 10 df
# using all the defaults
expsampleN.TOST(CV=0.3, dfCV=10)
# -> gives n=48 with achieved expected power 0.805082
# Compare this to the usual sample size with CV known as 'carved in stone'
sampleN.TOST(CV=0.3)

# More then one CV with corresponding degrees of freedom 
# other parameters as above
CVs <- c(0.25, 0.3)
dfs <- c(  22, 10)
expsampleN.TOST(CV=CVs, dfCV=dfs)
# -> gives a pooled CV=0.2664927 with df=32
# and a sample size n=34 with achieved expected power 0.815019
}