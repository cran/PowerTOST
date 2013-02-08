\name{sampleN.scABEL}
\alias{sampleN.scABEL}
\title{
Sample size estimation for BE via scaled (widened) BE acceptance limits
}
\description{
This function performs the Sample size estimation for the BE decision via 
scaled (widened) BE acceptance limits based on simulations.
}
\usage{
sampleN.scABEL(alpha = 0.05, targetpower = 0.8, theta0, theta1, theta2, 
               CV, design = c("2x3x3", "2x2x4"), regulator = c("EMA", "FDA"), 
               nsims = 1e+06, nstart, print = TRUE, details = TRUE)
}
\arguments{
  \item{alpha}{
Type I error probability. Per convention mostly set to 0.05.
}
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1.\cr
Typical values are 0.8 or 0.9.
}
  \item{theta0}{
'True' or assumed bioequivalence ratio. \cr
Defaults to 0.95 if given explicitly.
}
  \item{theta1}{
Lower limit for the point estimator constraint.\cr
Defaults to 0.8 if not given explicitly.
}
  \item{theta2}{
Upper limit for the point estimator constraint.\cr
Defaults to 1.25 if not given explicitly.
}
  \item{CV}{
Coefficient(s) of variation as ratio.\cr
If \code{length(CV) = 1} the same CV is assumed for Test and Reference.\cr
If \code{length(CV) = 2} the CV for Test must be given in CV[1] and for 
Reference in CV[2].
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
Number of simulations to be performed to obtain the (empirical) power.
}
  \item{nstart}{
Set this to a start for the sample size if a previous run failed.
}
  \item{print}{
If TRUE (default) the function prints its results.\cr 
If FALSE only the sample size will be returned. 
}
  \item{details}{
If set to \code{TRUE}, the default, the steps during sample size search are shown.
}
}
\details{
The simulations are done via the distributional properties of the statistical
quantities necessary for deciding BE based on widened ABEL.
}
\value{
Returns the numerical value of the sample size.
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
We are doing the sample size estimation only for balanced designs since the 
break down of the total subject number in case of unbalanced sequence groups
is not unique. Moreover the formulas used are only for balanced designs.\cr
In case of \code{regulator="FDA"} the sample size is only approximate since
the BE decision method is not exactly what is expected by the FDA. But the two
Laszlos state that the scABEL method should be 'operational' equivalent to the
FDA method. Thus the sample size should be comparable. 
}

\section{Warning }{
The sample size estimation for theta0 >= 1.2 and <=0.83 may be very time 
consuming and will eventually also fail since the start values chosen are
not really reasonable in that ranges. This is especially true if you choose
a target power >0.8 and your CV is >0.5.\cr
If you really need sample sizes in that range be prepared to restart the
sample size estimation via the argument nstart.
}

\seealso{
\code{\link{power.scABEL}}
}
\examples{
#using all the defaults
sampleN.scABEL(CV=0.3)
# results in a sample size n=27

# for the full replicate design, target power = 90%
# true assumed ratio = 0.9, FDA regulatory settings
sampleN.scABEL(CV=0.4, targetpower=0.9, theta0=0.9, design="2x2x4", regulator="FDA")
# results in a sample size n=32

}
