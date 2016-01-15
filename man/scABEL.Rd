\name{scABEL}
\alias{scABEL}
\title{
Helper function to calculate the (widened) scaled BE acceptance limits
}
\description{
The (widened) scaled BE acceptance limits are calculated according to the 
regulatory settings of EMA, FDA, ANVISA or via user defined constants.
}
\usage{
scABEL(CV, regulator = c("EMA", "ANVISA", "FDA", "USER"), r_const, CVswitch, CVcap)
}
\arguments{
  \item{CV}{
Coefficient of variation (of the Reference) as ratio
}
  \item{regulator}{
Regulatory body settings for the widening of the BE acceptance limits.\cr
Defaults to \code{regulator="EMA"}.
}
  \item{r_const}{
Regulatory constant for widening.
}
  \item{CVswitch}{
CV for switch to the widened ABEL.
}
  \item{CVcap}{
CV for cap of widening.
}
}
\details{
r_const, CVswitch and CVcap must be given if \code{regulator=="USER"}.\cr
Otherwise these arguments may be missing.
}
\value{
Returns a vector of lenghth 2 if one CV is given or a matrix if CV is given as vector 
with named components "lower" and "upper" of the scaled acceptance limits.
}
\author{
D. Labes
}
\seealso{
\code{\link{power.scABEL}}, \code{\link{sampleN.scABEL}}
}
\examples{
scABEL(CV=0.3, regulator="EMA")
# should give the usual limits:
# lower  upper
#  0.80   1.25
scABEL(CV=0.4, regulator="EMA")
# should give the widened limits:
#   lower     upper
# 0.746177 1.340165
}
