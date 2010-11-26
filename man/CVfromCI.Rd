\name{CVfromCI}
\alias{CVfromCI}
\title{CV from a given Confidence interval}
\description{
Calculates the CV (coefficient of variation) from a known confidence interval
of a BE study. Useful if no CV but CI was given in literature.
}
\usage{
CVfromCI(point, lower, upper, n, design = "2x2", alpha = 0.05)
}
\arguments{
  \item{point}{
Point estimator of the BE ratio.\cr
The point estimator can be missing.\cr
In that case it will be calculated as geometric mean of lower and upper.
}
  \item{lower}{
Lower confidence limit of the BE ratio.
}
  \item{upper}{
Upper confidence limit of the BE ratio.
}
  \item{n}{
Total number of subjects under study in case of a crossover study.\cr
Number of subjects per treatment in case of a parallel-group study.\cr
}
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
}
  \item{alpha}{
Error probability. Set it to (1-confidence)/2.\cr
Is 0.05 for the usual 90\% confidence intervals.
}
}
\details{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
}
\value{
Returns the value of the CV as ratio.
}
\author{D. Labes}
\note{
The calculations are based on the assumption of evaluation via log-transformed values.\cr
The calculations are further based on a common variance of Test and Reference \cr
treatments in replicate crossover studies or parallel group study, respectively.\cr
It is assumed that the sequence groups in a crossover study or the 
treatment arms in a parallel-group study are balanced.\cr
The estimated CV is conservative (i.e. greater than actually observed) in case of 
unbalanced studies.
}

\examples{
# Given a 90\% confidence interval (without point estimator) 
# from a classical 2x2 crossover with 22 subjects
CVfromCI(lower=0.91,upper=1.15,n=22,design="2x2")
# will give 
# [1] 0.2279405 i.e a CV ~ 23%
}