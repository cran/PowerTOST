\encoding{utf-8}
\name{power.TOST}
\alias{power.TOST}
\title{
Power of the classical TOST procedure
}
\description{
Calculates the exact or approximate power of the two-one-sided t-tests (TOST)
procedure for various study designs used in BE studies.
}
\usage{
power.TOST(alpha = 0.05, logscale = TRUE, theta1, theta2, theta0, CV, n, 
           design = "2x2", method="exact", robust=FALSE)
}
\arguments{
  \item{alpha}{
Significance level (one-sided). Commonly set to 0.05.
}
  \item{logscale}{
Should the data used on log-transformed or on original scale? \code{TRUE} (default) or \code{FALSE}.
}
  \item{theta0}{
\sQuote{True} or assumed T/R ratio or difference.\cr
In case of \code{logscale=TRUE} it must be given as ratio T/R.\cr
If \code{logscale=FALSE}, the difference in means. In this case, the difference may be expressed in two ways: relative to the same (underlying) reference mean, i.e. as (T-R)/R = T/R - 1; or as difference in means T-R. Note that in the former case the units of \code{CV}, \code{theta1} and \code{theta2} need also be given relative to the reference mean (specified as ratio).\cr
Defaults to 0.95 if \code{logscale=TRUE} or to 0.05 if \code{logscale=FALSE}
}
  \item{theta1}{
Lower (bio-)equivalence limit.\cr
In case of \code{logscale=TRUE} it is given as ratio.\cr
If \code{logscale=FALSE}, the limit may be expressed in two ways:
difference of means relative to the same (underlying) reference mean or in units of the difference of means.
Note that in the former case the units of \code{CV}, \code{theta0} and \code{theta2} need also be given relative to the reference mean (specified as ratio).\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
}
  \item{theta2}{
Upper (bio-)equivalence limit.\cr
In case of \code{logscale=TRUE} it is given as ratio.
If \code{logscale=FALSE}, the limit may be expressed in two ways:
difference of means relative to the same (underlying) reference mean or in units of the difference of means.
Note that in the former case the units of \code{CV}, \code{theta0} and \code{theta1} need also be given relative to the reference mean (specified as ratio).\cr
If not given, \code{theta2} will be calculated as \code{1/theta1} if \code{logscale=TRUE} or as \code{-theta1} if \code{logscale=FALSE}.
}
  \item{CV}{
In case of \code{logscale=TRUE} the (geometric) coefficient of variation given as ratio.\cr
If \code{logscale=FALSE} the argument refers to (residual) standard deviation of the response. In this case, standard deviation may be expressed two ways: relative to a reference mean (specified as ratio sigma/muR), i.e. again as a coefficient of variation; or untransformed, i.e. as standard deviation of the response. Note that in the former case the units of \code{theta0}, \code{theta1} and \code{theta2} need also be given relative to the reference mean (specified as ratio).\cr

In case of cross-over studies this is the within-subject CV, in case of a parallel-group design the CV of the total variability.
}
  \item{n}{
Number of subjects under study.\cr
Is total number if given as scalar, else number of subjects in the (sequence) 
groups. In the latter case the length of \code{n} vector has to be equal to the 
number of (sequence) groups.
}
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
  }
  \item{method}{
Method for calculation of the power.\cr  
Defaults to \code{"exact"} in which case the calculation is done based on formulas 
with \enc{Owen’s}{Owen's} Q. The calculation via \enc{Owen’s}{Owen's} Q can also be choosen with 
\code{method="owenq"}.\cr
Another exact method via direct integration of the bivariate non-central 
\emph{t}-distribution may be chosen with \code{method="mvt"}. This may have somewhat 
lower precision compared to \enc{Owen’s}{Owen's} Q and longer run-time.\cr
Approximate calculations can be choosen via \code{method="noncentral"} or 
\code{method="nct"} for the approximation using the non-central \emph{t}-distribution.
With \code{method="central"} or \code{method="shifted"} the relative crude
approximation via \sQuote{shifted} central \emph{t}-distribution is chosen.\cr
The strings for \code{method} may be abbreviated.
  }
  \item{robust}{
Defaults to \code{FALSE}. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the \sQuote{robust} 
evaluation (aka \enc{Senn’s}{Senn's} basic estimator). These degrees of freedom are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.\cr
Has only effect for higher-order crossover designs.
  }
}
\details{
The exact calculations of power are based on \enc{Owen’s}{Owen's} Q-function or by direct 
integration of the bivariate non-central \emph{t}-distribution via function
\code{\link[mvtnorm]{pmvt}} of package \code{mvtnorm}.\cr
Approximate power is implemented via the non-central \emph{t}-distribution
or the \sQuote{shifted} central \emph{t}-distribution.\cr\cr
The formulas cover balanced and unbalanced studies w.r.t (sequence) groups.\cr\cr
In case of parallel group design and higher order crossover designs 
(replicate crossover or crossover with more than two treatments) the calculations 
are based on the assumption of equal variances for Test and Reference products
under consideration.\cr\cr
The formulas for the paired means 'design' do not take a correlation parameter 
into account. They are solely based on the paired \emph{t}-test (TOST of differences = zero).
}
\value{
Value of power according to the input arguments.
}
\note{
Of course it is highly recommended to use the default \code{method="exact"} :-).\cr
There is no reason beside testing and for comparative purposes to use an 
approximation if the exact method is available.\cr\cr
If \code{n} is given as scalar (total sample size) and this number is not 
divisible by the number of (sequence) groups of the design an unbalanced design 
with small imbalance is assumed. A corresponding message is thrown showing the 
assumed numbers of subjects in (sequence) groups.\cr
The function does not vectorize properly if design is a vector.\cr
The function vectorizes properly if CV or theta0 are vectors.\cr
Other vector input is not tested yet.\cr\cr
The former function \code{power2.TOST()} designd to handle unbalanced studies is 
defunct since \code{power.TOST()} handles balanced as well as unbalanced designs.
}
\references{
Phillips KF. \emph{Power of the Two One-Sided Tests Procedure in Bioequivalence.} J Pharmacokin Biopharm. 1990;18(2):137--44. \doi{10.1007/BF01063556}

Diletti D, Hauschke D, Steinijans VW. \emph{Sample Size Determination for Bioequivalence Assessment by Means of 
Confidence Intervals.} Int J Clin Pharmacol Ther Toxicol. 1991;29(1):1--8. 
}
\author{
D. Labes, direct integration of bivariate non-central \emph{t}-distribution by B. Lang
}
\seealso{
 \code{\link{sampleN.TOST}, \link{known.designs}}
}
\examples{
# power for the 2x2 cross-over design with 24 subjects and CV 25\%
# using all the other default values
power.TOST(CV = 0.25, n = 24)
# should give: [1] 0.7391155
# nct approximation very good for this configuration
power.TOST(CV = 0.25, n = 24, method = "nct")
# gives also: [1] 0.7391155
# shifted-central-t  approximation 
power.TOST(CV = 0.25, n = 24, method = "shifted")
# gives:      [1] 0.7328894

# power for the 2x2 cross-over study with 24 subjects, CV 25\%
# and 2 drop-outs in the same sequence group (unbalanced study)
power.TOST(CV=0.25, n=c(10,12))
# should give: [1] 0.6912935
# not the same compared to the balanced setting
power.TOST(CV=0.25, n=22)
# should give: [1] 0.6953401
}
