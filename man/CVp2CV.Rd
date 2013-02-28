\name{CVp2CV}
\alias{CVp2CV}
\title{
Decompose CV(T) and CV(R) from 'pooled' CV of T/R 
}
\description{
Helper function to calculate CV(T) and CV(R) from a pooled CV(T/R)
assuming a ratio of the intra-subject variances. 
}
\usage{
CVp2CV(CV, ratio = 1.5)
}
\arguments{
  \item{CV}{
'pooled' CV of T/R.
}
  \item{ratio}{
Ratio of the intra-subject variances s^2(T)/s^2(R).
}
}
\details{
In case of knowing only the CV(T/R) f.i. from an ordinary cross-over you can 
calculate the components CV(T) and CV(R) assuming a ratio of the 
intra-subject variances.\cr
The formula the function is based on:\cr
\code{log(1.0 + CV^2) = (s2D + (sWT^2 + sWR^2)/2)}\cr
Insert sWT^2 = ratio* sWR^2 and solve for sWR^2.
The subject-by-formulation interaction term is assumed zero in this implementation.
}
\value{
Returns a numeric vector of the CV values for Test and Reference.
}
\author{
D. Labes
}
\examples{
CVp2CV(0.4, ratio=2)
# gives
# [1] 0.4677952 0.3225018
}
