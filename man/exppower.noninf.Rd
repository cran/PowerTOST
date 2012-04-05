\name{exppower.noninf}
\alias{exppower.noninf}
\title{
'Expected' power of non-inferiority test
}
\description{
Calculates the 'expected' power according to Julious for a variety of study
designs used in bioequivalence studies.
}
\usage{
exppower.noninf(alpha = 0.025, logscale=TRUE, theta0, margin,  
                CV, dfCV, n, design = "2x2", robust=FALSE)
}
\arguments{
  \item{alpha}{
Type I error probability, significance level. Defaults here to 0.025.
  }
  \item{logscale}{
Should the data used on log-transformed or on original scale? TRUE or FALSE.\cr
Defaults to TRUE.
  }  
  \item{theta0}{
'True' or assumed bioequivalence ratio or difference. 
Typically set to 0.95 (default if missing) if \code{logscale=TRUE}.\cr
Defaults to -0.05 if \code{logscale=FALSE}.
  }
  \item{margin}{
Non-inferiority margin.\cr
In case of \code{logscale=TRUE} it must be given as ratio, otherwise as diff.\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
  }
  \item{CV}{
Coefficient of variation as ratio. 
  }
  \item{dfCV}{
Degrees of freedom for the CV (error/residual degree of freedom).
}
  \item{n}{
Number of subjects to be planned (ntotal).
}
  \item{design}{
Character string describing the study design.
See \code{known.designs()} for designs covered in this package.
}
  \item{robust}{
Defaults to FALSE.\cr
Set to \code{TRUE} will use the degrees of freedom according to the 'robust' evaluation
(aka Senn's basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.
}
}
\details{
This function calculates the so-called 'expected' power based on formulas 
according to S.A. Julious. These take into account that usually the CV is not 
known but estimated from a previous study / studies with an uncertainty. 
See references.
}
\value{
Value of expected power according to the input.
}
\references{
S.A. Julious\cr
"Sample sizes for Clinical Trials"\cr
CRC Press, Chapman & Hall 2010
}
\author{
D. Labes
}

\seealso{
 \code{\link{expsampleN.noninf}, \link{power.noninf}, \link{power.TOST}} 
}
\examples{
# expected power for non-inferiority test of a 2x2 crossover
# CV 30\% known from a pilot study with 12 subjects (-> dfCV=10)
# using all the defaults for other parameters
# should give: [1] 0.6751358
exppower.noninf(CV=0.3, dfCV=10, n=40)

# Compare this to the usual power (CV known, "carved in stone")
# should give: [1] 0.7228685
power.noninf(CV=0.3, n=40)
}
