\name{PowerTOST-internal}
\title{Internal PowerTOST objects}
\alias{PowerTOSTinternal}
\alias{.approx.power.TOST}
\alias{.approx2.power.TOST}
\alias{.design.no}
\alias{.design.props}
\alias{.exppower.TOST}
\alias{.expsampleN0}
\alias{.power.TOST}
\alias{.power.RatioF}
\alias{.Q.integrand}
\alias{.sampleN0}
\alias{.sampleN0.RatioF}
\description{Internal PowerTOST objects: \cr
.approx.power.TOST: approximate power using non-central t-distribution\cr
.approx2.power.TOST: approximate power using central 'shifted' t-distribution\cr
.design.no: no of design from design literal\cr
.design.props: design properties (df, bk, step, name)\cr
.exppower.TOST: expected power\cr
.expsampleN0: start value of n for sample size based on exp. power.\cr
.power.TOST: power of classical TOST (working horse)\cr
.power.RatioF: power of ratio of means with normal data on original scale (working horse)\cr
.Q.integrand: Integrand for Owen's Q function \cr
.sampleN0: Start value of sample size\cr
.sampleN0.RatioF: Start value of sample size for ratio of means with normal data 
on original scale
}
\details{These are usually not to be called by the user.\cr
Due to the necessary NAMESPACE from R14.0 on these functions are no longer exported.}
\keyword{internal}