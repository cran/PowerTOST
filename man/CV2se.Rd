\name{CV2se/se2CV}
\alias{CV2se}
\alias{se2CV}
\title{Helper functions}
\description{Calculates the standard error from a given CV and 
vice versa for log-normal data}
\usage{CV2se(CV)
se2CV(se)}
\arguments{
  \item{CV}{coefficient of variation}
  \item{se}{standard error}
}
\value{
Returns se=sqrt(log(CV^2+1))\cr
or CV=sqrt(exp(se*se)-1)
}
\author{D. Labes}
\note{
These functions were intended for internal use.\cr
But may be usefull for others.
}

\examples{
CV2se(0.3)
#should give: [1] 0.2935604 

se2CV(0.2935604)
#[1] 0.3 
}

