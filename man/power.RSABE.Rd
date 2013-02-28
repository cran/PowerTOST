\name{power.RSABE}
\alias{power.RSABE}
\title{
(Empirical) Power for BE decision via linearized scaled ABE criterion
}
\description{
This function performs the power calculation of the BE decision via 
linearized scaled ABE criterion by simulations.
}
\usage{
power.RSABE(alpha = 0.05, theta1, theta2, theta0, CV, n, 
            design = c("2x3x3", "2x2x4"), regulator = c("FDA", "EMA"),
            nsims = 1e+05, details = FALSE, setseed=TRUE)
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. Conventionally mostly set to 0.05.
}
  \item{theta1}{
Conventional lower ABE limit to be applied in the mixed procedure if 
CVsWR <= CVswitch. Also lower limit for the point estimator constraint.\cr
Defaults to 0.8 if not given explicitly.
}
  \item{theta2}{
Conventional upper ABE limit to be applied in the mixed procedure if 
CVsWR <= CVswitch. Also upper limit for the point estimator constraint.\cr
Defaults to 1.25 if not given explicitly.
}
  \item{theta0}{
'True' or assumed bioequivalence ratio. \cr
Defaults to 0.95 if not given explicitly.
}
  \item{CV}{
Coefficient(s) of variation as ratio.\cr
If \code{length(CV) = 1} the same CV is assumed for Test and Reference.\cr
If \code{length(CV) = 2} the CV for Test must be given in CV[1] and for 
Reference in CV[2].
}
  \item{n}{
Number of subjects under study.\cr
May be given as vector. In that case it is assumed that n contains the number
of subjects in the sequence groups.\cr\cr
If n is given as single number (total sample size) and this number is not 
divisible by the number of sequences of the design an unbalanced design is assumed.
A corresponding message is thrown showing the numbers of subjects
}
  \item{design}{
Design of the study to be planned.\cr
2x3x3 is the partial replicate design (TRR/RTR/RRT).\cr
2x2x4 is the full replicate design with 2 sequences and 4 periods.\cr
Defaults to \code{design="2x3x3"}
}
  \item{regulator}{
Regulatory body settings for the scaled ABE criterion.\cr
Defaults to \code{design="FDA"}.\cr
Also the scaled ABE criterion is usually calculated with the FDA constant
r_const=log(1.25)/0.25 you can override this behavior to use the EMA setting
r_const=0.76 to avoid the discontinuity at CV=30\% and be more stringent.
}
  \item{nsims}{
Number of simulations to be performed to obtain the empirical power.
Defaults to 100 000 = 1e+5.
}
  \item{details}{
If set to \code{TRUE} the computational time is shown as well as the components
for the BE decision.\cr
p(BE-ABE) is the simulated probability for the conventional ABE test.
p(BE-SABEc) is the probability that the 95\% CI of the ABE criterion is <0.\cr
p(BE-PE) is the probability that the point estimate is within theta1 ... theta2.\cr
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo) random number 
generator. To avoid differences in power for different runs a set.seed() 
is issued if \code{setseed=TRUE}, the default. 
}
}
\details{
The linearized scaled ABE criterion is calculated according to the SAS code
given in the FDA progesterone guidance.\cr\cr
The simulations are done via the distributional properties of the statistical
quantities necessary for deciding BE based on scaled ABE criterion.\cr
The sample mean difference T-R is assumed normally distributed with 
mean=log(theta0) and sd=sqrt(bkni*sum(1/ni)*sqrt(s2D + (s2WT + s2WR)/2). 
s2WT and s2WR are obtained from the corresponding CV's via s2..=log(CV^2+1).
The subject-by-formulation interaction s2D is assumed in the current 
implementation =0. The ni are the number of subjects in the sequence groups
and bkni is the design constant in terms of ni.\cr
The sample variance (sd2s) of the difference T-R is obtained from the chi-squared 
distribution of sd2s*df/s2d with df degrees of freedom. s2d is the assumed
population (true) value.\cr
The sample variance (s2WRs) of the difference R-R is obtained from the chi-squared 
distribution of s2WRs*dfRR/s2WR with dfRR degrees of freedom. s2WR is the assumed
population (true) value of the within-subject variance of the Reference.
}
\value{
Returns the value of the empirical power.
}
\references{
FDA "Draft Guidance on Progesterone"\cr
Recommended Apr 2010; Revised Feb 2011\cr
http://www.fda.gov/downloads/Drugs/GuidanceComplianceRegulatoryInformation/Guidances/UCM209294.pdf

Laszlo Tothfalusi and Laszlo Endrenyi\cr
"Sample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs"\cr
J. Pharm. Pharmaceut. Sci. (www.cspsCanada.org) 15(1) 73 - 84, 2011

Tothfalusi L., Endrenyi L. and A. Garcia Arieta\cr
"Evaluation of Bioequivalence for Highly Variable Drugs with Scaled Average Bioequivalence"\cr
Clin. Pharmacokin. 48/11, 725-743 (2009)
}
\author{
D. Labes
}
\seealso{
\code{\link{sampleN.RSABE}}, \code{\link{power.scABEL}}
}
\examples{
# using all the defaults:
# design="2x3x3" -> partial replicate
# ABE limits, PE constraint 0.8-1.25
# true ratio =0.95, 1E+5 simulations
power.RSABE(CV=0.4, n=24)
# should give
# [1] 0.0.80864
#
# to explore the simulation error due to the state of the
# random number generator
power.RSABE(CV=0.4, n=24, setseed=FALSE)
# will give something like
# [1] 0.8081

}
