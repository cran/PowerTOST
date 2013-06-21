\name{sampleN.NTIDFDA}
\alias{sampleN.NTIDFDA}
\title{
Sample size estimation for BE decision via FDA method for narrow therapeutic index drugs (NTID's)
}
\description{
This function performs the Sample size estimation for the BE decision via 
FDA method for NTID's based on simulations. The study design is the full 
replicate design 2x2x4.
}
\usage{
sampleN.NTIDFDA(alpha = 0.05, targetpower = 0.8, theta0, theta1, theta2, CV, 
                nsims = 1e+05, nstart, print = TRUE, details = TRUE, 
                setseed = TRUE)
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
Attention! Defaults here to 0.975 if not given explicitly. The value was chosen 
nearer to 1 because the potency (contents) settings for NTID's are tightened 
by the FDA.
}
  \item{theta1}{
Conventional lower ABE limit to be applied in the FDA procedure.\cr
Defaults to 0.8 if not given explicitly.
}
  \item{theta2}{
Conventional upper ABE limit to be applied in the FDA procedure.\cr
Defaults to 1.25 if not given explicitly.
}
  \item{CV}{
Coefficient(s) of variation as ratio.\cr
If \code{length(CV) = 1} the same CV is assumed for Test and Reference.\cr
If \code{length(CV) = 2} the CV for Test must be given in CV[1] and for 
Reference in CV[2].
}
  \item{nsims}{
Number of simulations to be performed to obtain the empirical power.
Defaults to 100 000 = 1e+5.
}
  \item{nstart}{
Set this to a start value for the sample size if a previous run failed.\cr
May be missing.
}
  \item{print}{
If \code{TRUE} (default) the function prints its results.\cr 
If \code{FALSE} only the resulting dataframe will be returned. 
}
  \item{details}{
If set to \code{TRUE}, the default, the steps during sample size search are shown.
Moreover the details of the method settings are printed.
}
  \item{setseed}{
Simulations are dependent on the starting point of the (pseudo) random number 
generator. To avoid differences in power values for different runs a 
\code{set.seed(123456)} is issued if \code{setseed=TRUE}, the default. 
}
}
\details{
The linearized scaled ABE criterion is calculated according to the SAS code
given in the FDA Warfarine guidance. For deciding BE the study must pass that criterion, the conventional ABE test and additional the test that the ratio of 
sWT/sWR is < 2.5.\cr\cr
The simulations are done via the distributional properties of the statistical
quantities necessary for deciding BE based on these method.\cr
Details can be found in a document "Implementation_scaledABE_sims" located in 
the doc subdirectory of the package.
}
\value{
Returns now a data.frame with the input and sample size results.\cr
The "Sample size" column contains the total sample size.\cr
The "nlast" column contains the last n value. May be useful for re-starting.
}
\references{
FDA "Draft Guidance on Warfarin Sodium"\cr
Recommended Dec 2012\cr
http://www.fda.gov/downloads/Drugs/GuidanceComplianceRegulatoryInformation/Guidances/UCM201283.pdf}
\author{
D. Labes
}
\note{
The design recommended by the FDA is the full replicate design 2x2x4. 
Only this design is implemented.\cr
The sample size estimation is done only for balanced studies since the 
break down of the total subject number in case of unbalanced sequence groups
is not unique. Moreover the formulas used are only valid for balanced designs.\cr
The FDA method is described for the ABE limits 0.8 ... 1.25 only. Setting theta1, theta2
to other values may not be reasonable and is not tested.
}
\section{Warning }{
For some input constellations the sample size search may be very time 
consuming and will eventually also fail since the start values chosen may
not really reasonable for them. This applies especially for theta0 values
near to the implied scaled (tightened/widened) ABE limits according to
\code{exp(+-1.053605*swR)}.\cr
In case of a failed sample size search you may restart with setting tha argument \code{nstart}.\cr
In case of theta0 values outside the implied scaled (tightened/widened) ABE limits
no sample size estimation is possible and the function throws an error 
(f.i. CV=0.04, theta0=0.95).
}

\seealso{
\code{\link{power.NTIDFDA}}
}
\examples{
sampleN.NTIDFDA(CV=0.04,theta0=0.975)
# should give
# n=54 with an (empirical) power of 0.809590
#
# Test formulation with lower variability
sampleN.NTIDFDA(CV=c(0.04,0.06),theta0=0.975)
# should give
# n=20 with an (empirical) power of 0.0.814610
}

