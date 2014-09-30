\name{sampleN.dp}
\alias{sampleN.dp}
\title{
Sample size estimation of dose-proportionality studies evaluated via Power model
}
\description{
Performes a sample size estimation for dose-proportionality studies using the Power model for 
cossover (Latin square) or parallel group designs via a confidence interval
equivalence criterion.
}
\usage{
sampleN.dp(alpha = 0.05, CV, doses, targetpower = 0.8, beta0 = 1, theta1 = 0.8, 
           theta2 = 1/theta1, design = c("crossover", "parallel"), 
           print = TRUE, details = FALSE, imax = 100)
}
\arguments{
  \item{alpha}{
Type 1 error. Usually taken as 0.05.
}
  \item{CV}{
Coefficient of variation. Is intra-subject CV for \code{design="crossover"} and
CV of total variability in case of \code{design="parallel"}
}
  \item{doses}{
Vector of dose values under study. At least 2 doses have to be given.
}
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1.\cr
Typical values are 0.8 or 0.9.
}
  \item{beta0}{
'True' slope of power model. Defaults to 1 also a real planing should be done 
with some deviation to 1, but within slope acceptance range according to
\code{1+log(theta1)/log(rd)} and \code{1+log(theta2)/log(rd)}
were rd is the ratio of highest to lowest dose. Otherwise the function issues 
an error.
}
  \item{theta1}{
Lower acceptance limit for the ratio of dose normalized means (Rdmn).\cr
Transformes into slope acceptance range as described under item \code{beta0}.
}
  \item{theta2}{
Upper acceptance limit for the ratio of dose normalized means (Rdmn).
}
  \item{design}{
Crossover design (default) or parallel group design.\cr
Crossover design means Latin square design with number of doses as dimension.
}
  \item{print}{
If \code{TRUE} (default) the function prints its results.\cr 
If \code{FALSE} only the data.frame with the results will be returned. 
}
  \item{details}{
If \code{TRUE} the design characteristics and the steps during
sample size calculations will be shown. \cr 
Defaults to \code{FALSE}.
}
  \item{imax}{
Maximum number of steps in sample size search. \cr 
Defaults to 100. Adaption only in rare cases if any needed.
}
}
\details{
The sample size is calculated via iterative evaluation of \code{power.dp()}.\cr
Start value for the sample size search is taken from a large sample approximation.\cr
The sample size is bound to 2*number of doses as minimum.\cr
Balanced designs are used although this is not absolutely necessary.
}
\value{
A data.frame with the input and results will be returned.\cr
The "Sample size" column contains the total sample size.
}
\references{
Patterson, Jones\cr
"Bioequivalence and Statistics in Clinical Pharmacology"\cr
Chapman & Hall/CRC, Boca Raton, 2006, page 239\cr

Hummel J, McKendrick S, Brindley C, and R French\cr
"Exploratory assessment of dose proportionality: review of current approaches\cr 
and proposal for a practical criterion"\cr
Pharmaceut Statist 8(1), 38-49 (2009)

Sethuraman VS, Leonov S, Squassante L, Mitchell TR, Hale MD\cr
"Sample size calculation for the Power Model for dose proportionality studies"\cr
Pharm Stat. Vol. 6(1):35-41 (2007)
}
\author{
D. Labes
}
\section{Warning }{
This function is 'experimental' only since there is a pending discussion of 
differences to the SAS code given in Patterson, Jones.\cr
Additionally the function is not thorougly tested yet.
}
\seealso{
\code{\link{power.dp}}
}
\examples{
# using all the defaults, i.e. crossover design, alpha=0.05
# theta1=0.8, theta2=1.25
sampleN.dp(CV=0.2, doses=c(1, 2, 8), beta0=1.02)
# should give n=18
}
