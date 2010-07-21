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
definit integral in Owen's Q.
}
\value{
Returns the value of Owen's Q.
}
\references{
Owen, D. B. (1965) \cr
"A Special Case of a Bivariate Non-central t-Distribution"\cr
Biometrika, 52, 437-446.
}
\author{D. Labes}
\note{
This function is intended for internal use in the power calculations.\cr
But may be usefull for others.
}

\examples{
#This function is intended for internal use.
OwensQ(10,2.5,5,0,2)
#[1] 9.388137e-06 
}

