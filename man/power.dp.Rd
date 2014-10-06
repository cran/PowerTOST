\name{power.dp}
\alias{power.dp}
\title{
Power of dose-proportionality studies evaluated via Power model
}
\description{
Calculates the power of dose-proportionality studies using the Power model for 
cossover (Latin square) or parallel group designs via a confidence interval
equivalence criterion.
}
\usage{
power.dp(alpha = 0.05, CV, doses, n, beta0 = 1, theta1 = 0.8, theta2 = 1/theta1, 
         design = c("crossover", "parallel"))
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
Vector of dose values. At least 2 doses have to be given.
}
  \item{n}{
Number of subjects. Is total number if given as scalar, else number of subjects
in the (sequence) groups. In the latter case the length of n vector has to be 
the same as length of vector doses.
}
  \item{beta0}{
'True' slope of power model.
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
}
\details{
The power calculations are based on TOST for testing equivalence of the slope = 1.\cr
Power is calculated via non-central t approximation only.
}
\value{
Value of power according to the input arguments.
}
\references{
Patterson S, Jones B\cr
"Bioequivalence and Statistics in Clinical Pharmacology"\cr
Chapman & Hall/CRC, Boca Raton, 2006, page 239\cr
(contains presumably an error)

Sethuraman VS, Leonov S, Squassante L, Mitchell TR, Hale MD\cr
"Sample size calculation for the Power Model for dose proportionality studies"\cr
Pharm. Stat. 2007 Jan-Mar;6(1):35-41

Hummel J, McKendrick S, Brindley C, and R French\cr
"Exploratory assessment of dose proportionality: review of current approaches\cr 
and proposal for a practical criterion"\cr
Pharm. Stat. 2009 Jan-Mar;8(1):38-49
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
\code{\link{sampleN.dp}}
}
\examples{
# using all the defaults, i.e. crossover design, alpha=0.05, beta0=1
# theta1=0.8, theta2=1.25
power.dp(CV=0.2, doses=c(1, 2, 8), n=12)
}
