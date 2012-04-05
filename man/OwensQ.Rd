\name{OwensQ}
\alias{OwensQ}
\title{Owen's Q-function}
\description{Calculates Owen's Q function.}
\usage{OwensQ(nu, t, delta, a, b)}
\arguments{
  \item{nu}{degree of Owen's Q}
  \item{t}{parameter t}
  \item{delta}{parameter delta}
  \item{a}{lower integration limit}
  \item{b}{upper integration limit}
}
\details{
Uses integrate() from package stats to perform the numerical evaluation of the  
definite integral in Owen's Q. \cr
See \url{../doc/BE_power_sample_size_excerpt.pdf} in the package sub-directory /doc 
for the definition of Owen's Q. \cr
The arguments to the function must be scalars. No vectors allowed.
}
\value{
Returns the value of Owen's Q-function at given input arguments.
}
\references{
Owen, D. B. (1965) \cr
"A Special Case of a Bivariate Non-central t-Distribution"\cr
Biometrika, 52, 437-446.
}
\author{D. Labes}
\note{
This function is intended for internal use in the power calculations.\cr
But may be useful for others.
}

\examples{
# This function is intended for internal use.
OwensQ(10,2.5,5,0,2)
#should give [1] 9.388137e-06 
}

