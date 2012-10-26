\name{CVCL}
\alias{CVCL}
\title{
Confidence limits of a CV for log-normal data
}
\description{
The function calculates the 1-alpha confidence limits (either 1-sided
or 2-sided) via the chi-squared distribution of the error variance the CV 
is based on.
}
\usage{
CVCL(CV, df, side = c("upper", "lower", "2-sided"), alpha = 0.05)
}
\arguments{
  \item{CV}{Coefficient of variation}
  \item{df}{degrees of freedom of the CV (error variance)}
  \item{side}{Side(s) to calculate the confidence limits for}
  \item{alpha}{Type I error probability, aka significance level}
}
\value{
Numeric value or numeric vector of the confidence limits.
}
\author{D. Labes}
\examples{
# upper one-sided 95\% CL of a CV=0.3 
# from a study with df=22 (f.i. a 2x2 crossover with n=24)
# side="upper" is standard if not explicitly given
CVCL(0.3, df=22)
# should give 0.4075525
}
