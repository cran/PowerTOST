\name{expsampleN.TOST}
\alias{expsampleN.TOST}
\title{
Sample size based on 'expected' power
}
\description{
Calculates the sample size based on Julious 'expected' power for a
variety of study designs used in bioequivalence studies. See known.designs() 
for the study designs covered
}
\usage{
expsampleN.TOST(alpha = 0.05, targetpower = 0.8, theta0 = 0.95, theta1 = 0.8, 
                theta2, CV, dfCV, alpha2 = 0.05, design = "2x2", robust=FALSE,
                print = TRUE, details = FALSE, imax=100)
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
  \item{robust}{
Defaults to FALSE. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' evaluation
(aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.
}
  \item{print}{
If TRUE (default) the function prints its results.\cr
If FALSE only a data.frame with the results will be returned.
}
  \item{details}{
If TRUE the design characteristics and the steps during
sample size calculations will be shown.\cr Defaults to FALSE.
}
  \item{imax}{
Maximum number of steps in sample size search.\cr 
Defaults to 100. Adaption only in very rare cases needed.\cr
Never seen a need for adaption up to now.
}
}
\details{
The sample size is calculated based on 'expected' power via Julious formulas 
based on non-central t-distribution.
}
\value{
A data.frame with the input values and the result of the sample size estimation.\cr
The "Sample size" column contains the \bold{total} sample size in case of the 
crossover designs and the paired means design.\cr
It contains also the \bold{total} sample size  in case of the 2-parallel groups design
from package version >0.9-0 on.
}
\references{
S.A. Julious, R.J. Owen\cr
"Sample size calculations for clinical studies allowing for 
uncertainty in variance"\cr
Pharmaceutical Statistics (2006), 5, 29-37

S.A. Julious\cr
"Sample sizes for Clinical Trials"\cr
CRC Press, Chapman & Hall, Boca Raton 2010

S. Senn\cr
"Cross-over Trials in Clinical Research"
Second edition\cr
Wiley, Chichester 2002
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