\name{power.scABEL}
\alias{power.scABEL}
\title{
(Empirical) Power for BE via scaled (widened) BE acceptance limits
}
\description{
This function performs the power calculation of the BE decision via 
scaled (widened) BE acceptance limits by simulations.
}
\usage{
power.scABEL(alpha = 0.05, theta1, theta2, theta0, CV, n, 
             design = c("2x3x3", "2x2x4"), regulator = c("EMA", "FDA"), 
             nsims = 1e+06, details = FALSE)
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. Conventionally mostly set to 0.05.
}
  \item{theta1}{
Lower limit for the point estimator constraint.\cr
Defaults to 0.8 if not given explicitly.
}
  \item{theta2}{
Upper limit for the point estimator constraint.\cr
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
Regulatory body settings for the widening of the BE acceptance limits.\cr
Defaults to \code{design="EMA"}
}
  \item{nsims}{
Number of simulations to be performed to obtain the empirical power.
}
  \item{details}{
If set to \code{TRUE} the computational time is shown as well as the components
for the BE decision.\cr
p(BE-wABEL) is the probability that the CI is within widened limits.\cr
p(BE-PE) is the probability that the point estimate is within theta1 ... theta2.\cr
p(BE-ABE) is the simulated probability for the conventional ABE test.
 
}
}
\details{
The simulations are done via the distributional properties of the statistical
quantities necessary for deciding BE based on widened ABEL.
}
\value{
Returns the value of the empirical power.
}
\references{
Laszlo Tothfalusi and Laszlo Endrenyi
"Sample Sizes for Designing Bioequivalence Studies for Highly Variable Drugs"
J. Pharm. Pharmaceut. Sci. (www.cspsCanada.org) 15(1) 73 - 84, 2011
}
\author{
D. Labes
}
\note{
In case of \code{regulator="FDA"} the (empirical) power is only approximate since
the BE decision method is not exactly what is expected by the FDA. But the two
Laszlos state that the scABEL method should be 'operational' equivalent to the
FDA method.
}

\seealso{
\code{\link{sampleN.scABEL}}
}
\examples{
# using all the defaults:
# design="2x3x3", EMA regulatory settings
# PE constraint 0.8-1.25, cap on widening if CV>0.5
# true ratio =0.95, 1E+6 simulations
power.scABEL(CV=0.4, n=29)
# Unbalanced design. n(i)=10/10/9 assumed.
# [1] 0.81037
# with details=TRUE to view the computational time
power.scABEL(CV=0.5, n=54, theta0=1.15, details=TRUE)
#should give
# 1e+06 sims. Time elapsed (sec):
#   user  system elapsed 
#   0.80    0.07    0.86 
# p(BE-ABE)= 0.274807 ; p(BE-wABEL)= 0.817483 ; p(BE-PE)= 0.855218 
#[1] 0.807211
}
