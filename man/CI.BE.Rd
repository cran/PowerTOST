\name{CI.BE}
\alias{CI.BE}
\title{
1-2*alpha confidence interval given point estimator, CV and n
}
\description{
Utility function to calculate the 1-2*alpha CI's given point est., CV and n for
the various designs covered in this package.
}
\usage{
CI.BE(alpha = 0.05, pe, CV, n, design = "2x2", robust = FALSE)
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. Defaults to 0.05.
}
  \item{pe}{
Point estimator (GMR).
}
  \item{CV}{
Coefficient of variation of error variability as ratio. 
}
  \item{n}{
Total number of subjects if a scalar is given.\cr
Number of subjects in (sequence) groups if given as vector.
}
  \item{design}{
Character string describing the study design.\cr
See \code{known.designs()} for designs covered in this package.
}
  \item{robust}{
Defaults to FALSE.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' 
evaluation (aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.
}
}
\value{
Returns the 1-2*alpha confidence interval.\cr
Returns a vector with named elements lower, upper if arguments pe and CV are 
scalars, else a matrix with columns lower, upper is returned.
}
\author{
D. Labes
}
\note{
The function assumes an evaluation using log-transformation.\cr
The function assumes equal varainces in case of \code{design="parallel"} and 
the higher order crossover designs.\cr
The formula implemented covers balanced and unbalanced designs.\cr\cr
If the function vectorizes propperly is not thoroughly tested.
}

\examples{
# 90\% confidence interval for the 2x2 crossover
# n(total) = 24
CI.BE(pe=0.95, CV=0.3, n=24)
# should give
#    lower     upper 
#0.8213465 1.0988055 
# same with number of subjects in sequence groups
CI.BE(pe=0.95, CV=0.3, n=c(12, 12))
}