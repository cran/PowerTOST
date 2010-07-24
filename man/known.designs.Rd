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
This function is for informal purposes and will be used internal
for obtaining characteristics if the designs used in calculation formulas.
}
\value{
Returns a data.frame with \cr
  no = number of the design, \cr
  design = string for identifying the design, \cr
  df = degrees of freedom of the design,\cr
  steps = step width in the iterative sample size estimation, \cr
  bk = so-called design constant (here in terms of total n)
}
\note{
The df for the replicate designs are those without carry-over in the model.\cr
Chen, Chow and Liu used models with carry-over, i.e. one df lower as here.\cr
n is the total number of subjects in case of cross-over design or 
the number per group in case of parallel group design.\cr
df2 = degrees of freedom for robust analysis (aka Senn's basic estimator)
are not used at moment.
}
\references{
K.-W. Chen, S.-C. Chow and G. Liu\cr
"A Note on Sample Size Determination for Bioequivalence Studies\cr
with Higher-order Crossover Designs" \cr
J. Pharmacokinetics and Biopharmaceutics, Vol. 25, No. 6, p753-765 (1997)

S. Senn\cr
"Cross-over Trials in Clinical Research" \cr
Second Edition, John Wiley & Sons, Chicester 2002
}
\author{
D.Labes
}
\examples{
known.designs()
}
