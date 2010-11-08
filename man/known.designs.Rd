\name{known.designs}
\alias{known.designs}
\title{Show the 'known' designs}
\description{
Returns the known study designs for which power and sample size
can be calculated within this package.
}
\usage{
known.designs()
}
\details{
This function is for informal purposes and will be used internal \cr
for obtaining characteristics if the designs used in calculation formulas.
}
\value{
Returns a data.frame with \cr
  no = number of the design, \cr
  design = string for identifying the design, \cr
  df = degrees of freedom of the design,\cr
  steps = step width in the iterative sample size estimation, \cr
  bk = so-called design constant (here in terms of total n)\cr\cr
The design string has to be used in the functions call for power and sample size.
}
\note{
The design string for higher order crossover designs is named as:\cr
\code{treatments x sequences x periods} in case of replicate designs and\cr
\code{treatments x periods} in case of crossover designs for more then 2 treatments.\cr
\cr
The df for the replicate crossover designs are those without carry-over in the model.\cr
Chen, Chow and Liu used models with carry-over, i.e. one df lower as here.\cr
The design constant in case of design 2x2x4 is here bk=1. 
Chen, Chow and Liu used bk=1.1 due to carry-over in the model.\cr\cr
n is the total number of subjects in case of cross-over design or 
the number per treatment group in case of parallel group design.\cr\cr
df2 = degrees of freedom for robust analysis (aka Senn's basic estimator). 
These are not used at moment.
}
\references{
K.-W. Chen, S.-C. Chow and G. Liu\cr
"A Note on Sample Size Determination for Bioequivalence Studies with \cr
Higher-order Crossover Designs" \cr
J. Pharmacokinetics and Biopharmaceutics, Vol. 25, No. 6, p753-765 (1997)

S. Senn\cr
"Cross-over Trials in Clinical Research" \cr
Second Edition, John Wiley & Sons, Chichester 2002
}
\author{
D. Labes
}
\examples{
known.designs()
}
