\name{CVfromCI}
\alias{CVfromCI}
\title{CV from a given Confidence interval}
\description{
Calculates the CV (coefficient of variation) from a known confidence interval
of a BE study.\cr 
Useful if no CV but the 90\% CI was given in literature.
}
\usage{
CVfromCI(point, lower, upper, n, design = "2x2", alpha = 0.05, robust=FALSE)
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
Total number of subjects under study.
  }
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
  }
  \item{alpha}{
Error probability. Set it to \code{(1-confidence)/2}.\cr
Is 0.05 for the usual 90\% confidence intervals.
  }
  \item{robust}{
With \code{robust=FALSE} the usual degrees of freedom of the designs are used.\cr
With \code{robust=TRUE} the degrees of freedom for the so-called robust
evaluation (df2 in known.designs()) will be used. This may be helpful
if the CI was evaluated via mixed model or via intra-subject contrasts 
(aka Senn's basic estimator).
  }
}
\details{
See Helmut Schuetz lectures at \url{www.bebac.at/lectures.htm} for a description
of the algebra underlying this function.
}
\value{
Numeric value of the CV as ratio.
}
\author{D. Labes}
\note{
The calculations are based on the assumption of evaluation via log-transformed values.\cr
The calculations are further based on a common variance of Test and Reference 
treatments in replicate crossover studies or parallel group study, respectively.\cr
It is assumed that the sequence groups in a crossover study or the 
treatment arms in a parallel-group study are balanced.\cr
The estimated CV is conservative (i.e. greater than actually observed) in case of 
unbalanced studies.
}

\examples{
# Given a 90\% confidence interval (without point estimator) 
# from a classical 2x2 crossover with 22 subjects
CVfromCI(lower=0.91, upper=1.15, n=22, design="2x2")
# will give 
# [1] 0.2279405 i.e a CV ~ 23%
}