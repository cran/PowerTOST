\name{power.2TOST}
\alias{power.2TOST}
\title{
Power for two simultaneous TOST procedures
}
\description{
Calculates the exact power of two simultaneous TOST procedures (where
the two parameters of the two TOSTs are correlated with some correlation)
for various study designs used in BE studies 
}
\usage{
power.2TOST(alpha = c(0.05, 0.05), logscale = TRUE, theta1, theta2, theta0, 
            CV, n, rho, design = "2x2", robust = FALSE, setseed = TRUE)
}
\arguments{
  \item{alpha}{
  Vector; contains one-sided significance level for each of the two TOSTs.\cr
  For one TOST, by convention mostly set to 0.05.
}
  \item{logscale}{
Should the data used on log-transformed or on original scale? TRUE or FALSE.\cr
Defaults to TRUE.
}
  \item{theta1}{
Vector; contains lower bioequivalence limit for each of the two TOSTs.\cr
In case of \code{logscale=TRUE} it is given as ratio, otherwise as diff. to 1.\cr
Defaults to \code{c(0.8, 0.8)} if \code{logscale=TRUE} or to \code{c(-0.2, -0.2)}
if \code{logscale=FALSE}.
}
  \item{theta2}{
Vector; contains upper bioequivalence limit for each of the two TOSTS.\cr
If not given theta2 will be calculated as \code{1/theta1} if \code{logscale=TRUE} \cr
or as \code{-theta1} if \code{logscale=FALSE}.
}
  \item{theta0}{
Vector; contains 'true' assumed bioequivalence ratio for each of the two TOSTs. \cr
In case of \code{logscale=TRUE} each element must be given as ratio, \cr
otherwise as difference to 1. See examples.\cr
Defaults to \code{c(0.95, 0.95)} if \code{logscale=TRUE} or to 
\code{c(0.05, 0.05)} if \code{logscale=FALSE}.
}
  \item{CV}{
Vector of coefficient of variations (given as as ratio, e.g. 0.2 for 20\%).\cr
In case of cross-over studies this is the within-subject CV, \cr
in case of a parallel-group design the CV of the total variability.\cr
In case of \code{logscale=FALSE} CV is assumed to be the respective standard 
deviation.
}
  \item{n}{
Number of subjects under study.\cr
Is total number if given as scalar, else number of subjects in the (sequence) 
groups. In the latter case the length of n vector has to be equal to the 
number of (sequence) groups.
}
  \item{rho}{
Correlation between the two parameters under consideration. This is defined 
as correlation between the estimator of the treatment difference of
parameter one and the estimator of the treatment difference of parameter two.
}
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
}
  \item{robust}{
Defaults to FALSE. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' 
evaluation (aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.\cr
Has only effect for higher-order crossover designs.
  }
  \item{setseed}{
Calculation depends on \code{pmvt()} which is based on randomized quasi Monte Carlo
methods. If \code{setseed=TRUE} a seed value is set, the default.
}
}
\details{
The exact calculations of the power are performed via integration of the
4-dimensional non-central t-distribution via function \code{pmvt()} of package 
\code{mvtnorm}. An absolute error tolerance of 1e-04 is set within 
\code{pmvt()}.\cr\cr
The formulas cover balanced and unbalanced studies w.r.t (sequence) groups.\cr\cr
In case of parallel group design and higher order crossover designs 
(replicate crossover or crossover with more than two treatments) the calculations 
are based on the assumption of equal variances for Test and Reference products
under consideration.\cr\cr
The formulas for the paired means 'design' do not take an additional correlation
parameter into account. They are solely based on the paired t-test 
(TOST of differences = zero).
}
\value{
Value of power.
}
\note{
If n is given as scalar (total sample size) and this number is not 
divisible by the number of (sequence) groups of the design an unbalanced design 
with small imbalance is assumed. A corresponding message is thrown showing the 
assumed numbers of subjects in (sequence) groups.\cr
The function does not vectorize propperly if design is a vector. Moreover,
theta0 and CV must be of length two, thus further vectorizing is not possible.\cr
Other vector input is not tested yet.\cr\cr
}
\references{
Hua S. Y., Xu S., and D'Agostino Sr. R. B.\cr
"Multiplicity adjustments in testing for bioequivalence"\cr
Statistics in Medicine, Vol. 34, Issue 2, 215-231 (2015)

Lang B., Fleischer F. (in press).\cr
"Letter to the Editor 'Comments on Multiplicity adjustments in testing for bioequivalence'".\cr
Statistics in Medicine.

Phillips, K. F. (2009) \cr 
"Power for Testing Multiple Instances of the Two One-Sided Tests Procedure"\cr
The International Journal of Biostatistics: Vol. 5: Iss. 1, Article 15.
}
\author{
Benjamin Lang\cr
}
\seealso{
 \code{\link{sampleN.2TOST}, \link{known.designs}}
}
\examples{
# Power for the 2x2x2 cross-over design with 24 subjects, intra-subject
# standard deviation of 0.3 (CV = 30.7\%) and assumed ratios of 1.05 for both
# parameters, and correlation 0.75 between parameters (using all the other
# default values)
power.2TOST(theta0 = rep(1.05, 2), CV = rep(se2CV(0.3), 2), n = 24, rho = 0.75)
# should give: [1] 0.3981441

# Setting as before but use rho = 1 to replicate result of power.TOST()
p1 <- power.2TOST(theta0 = rep(1.05, 2), CV = rep(se2CV(0.3), 2), n = 24, rho = 1)
p2 <- power.TOST(theta0 = 1.05, CV = se2CV(0.3), n = 24)
all.equal(p1, p2, tolerance = 1e-04)
}