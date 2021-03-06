\encoding{utf-8}
\name{expsampleN.TOST}
\alias{expsampleN.TOST}
\title{
Sample size based on expected power
}
\description{
Estimates the sample size based on the expected power for a variety of study 
designs used in bioequivalence studies. See \link{known.designs} for the study 
designs covered.
}
\usage{
expsampleN.TOST(alpha = 0.05, targetpower = 0.8, logscale=TRUE, theta0, 
                theta1, theta2, CV, design = "2x2", robust = FALSE, 
                prior.type = c("CV", "theta0", "both"), prior.parm = list(),
                method = c("exact", "approx"), print = TRUE, details)
}
\arguments{
  \item{alpha}{
Significance level (one-sided). Commonly set to 0.05.
  }
  \item{targetpower}{
Power to achieve at least. Must be >0 and <1.
Typical values are 0.8 or 0.9.
  }
  \item{logscale}{
Should the data used on log-transformed or on original scale? \code{TRUE} (default) or \code{FALSE}.
  }  
  \item{theta0}{
Assumed \sQuote{true} (or \sQuote{observed} in case of \code{prior.type != "CV"}) bioequivalence 
ratio or difference.\cr
In case of \code{logscale=TRUE} it must be given as ratio T/R.\cr
If \code{logscale=FALSE}, the difference in means. In this case, the difference may be expressed in two ways: relative to the same (underlying) reference mean, i.e. as (T-R)/R = T/R - 1; or as difference in means T-R. Note that in the former case the units of \code{CV}, \code{theta1} and \code{theta2} need also be given relative to the reference mean (specified as ratio).\cr
Defaults to 0.95 if \code{logscale=TRUE} or to 0.05 if \code{logscale=FALSE}
  }
  \item{theta1}{
Lower (bio-)equivalence limit.\cr
In case of \code{logscale=TRUE} it is given as ratio.\cr
If \code{logscale=FALSE}, the limit may be expressed in two ways:
difference of means relative to the same (underlying) reference mean or in units of the difference of means.
Note that in the former case the units of \code{CV}, \code{theta0} and \code{theta2} need also be given relative to the reference mean (specified as ratio).\cr
Defaults to 0.8 if \code{logscale=TRUE} or to -0.2 if \code{logscale=FALSE}.
  }
  \item{theta2}{
Upper (bio-)equivalence limit.\cr
In case of \code{logscale=TRUE} it is given as ratio.
If \code{logscale=FALSE}, the limit may be expressed in two ways:
difference of means relative to the same (underlying) reference mean or in units of the difference of means.
Note that in the former case the units of \code{CV}, \code{theta0} and \code{theta1} need also be given relative to the reference mean (specified as ratio).\cr
If not given, \code{theta2} will be calculated as \code{1/theta1} if \code{logscale=TRUE} or as \code{-theta1} if \code{logscale=FALSE}.
  }
  \item{CV}{
In case of \code{logscale=TRUE} the (geometric) coefficient of variation given as ratio.\cr
If \code{logscale=FALSE} the argument refers to (residual) standard deviation of the response. In this case, standard deviation may be expressed two ways: relative to a reference mean (specified as ratio sigma/muR), i.e. again as a coefficient of variation; or untransformed, i.e. as standard deviation of the response. Note that in the former case the units of \code{theta0}, \code{theta1} and \code{theta2} need also be given relative to the reference mean (specified as ratio).\cr

If \code{prior.type="CV"} may be given as vector: The CVs are then pooled (as a 
weighted mean with their degrees of freedoms as weights).\cr

In case of cross-over studies this is the within-subject CV, in case of a parallel-group design the CV of the total variability.
  }
  \item{design}{
Character string describing the study design.\cr
See \link{known.designs} for designs covered in this package.
  }
  \item{robust}{
Defaults to FALSE. With that value the usual degrees of freedom will be used.\cr
Set to \code{TRUE} will use the degrees of freedom according to the \sQuote{robust} evaluation
(aka Senn\enc{’}{'}s basic estimator). These df are calculated as \code{n-seq}.\cr
See \code{known.designs()$df2} for designs covered in this package.
  }
  \item{prior.type}{
Specifies which parameter uncertainty should be accounted for. In case of 
\code{prior.type = "CV"} (the default), only the uncertainty with respect to the
CV will be considered (\emph{i.e.}, the given treatment effect is assumed to be fix).
In case of \code{prior.type = "theta0" } only uncertainty with respect to the
treatment ratio/difference will be accounted for (\emph{i.e.}, the given CV is assumed
to be fix). In case of  \code{prior.type = "both"} the power value will be
unconditional with respect to both the \code{CV} and \code{theta0}.
  }
  \item{prior.parm}{
A list of parameters expressing the prior information about the 
variability and/or treatment effect. Possible components are \code{df}, 
\code{SEM}, \code{m}, \code{design}.\cr
For \code{prior.type = "CV"} the degrees of freedom from the prior trial are
required. This information can be provided by specifying the single component 
\code{df} or the combination consisting of \code{m} and \code{design}.\cr
For \code{prior.type = "theta0"} the standard error of the treatment difference 
from the prior trial is required. This information can be provided by specifying
the single component \code{SEM} or the combination consisting of 
\code{m} and \code{design}.\cr
For \code{prior.type = "both"} the degrees of freedom and the standard error of
the treatment difference are required. This information can be provided by 
specifying the combination consisting of \code{df} and \code{SEM} or via the 
combination \code{m} and \code{design}.\cr
See 'Details' for a technical description on each component.
  }
  \item{method}{
Defaults to \code{method="exact"}.
In that case the expected power will be calculated as expected value of the
power with respect to the (prior) distribution of the respective parameter(s).\cr
Set to \code{method="approx"} the expected power according to the 
approximate formulas given in the book from Julious or in the Julious/Owen 
paper will be calculated (using non-central \emph{t}); this only affects 
\code{prior.type = "CV"}.
}
  \item{print}{
If \code{TRUE} (default) the function prints its results. If \code{FALSE} only a data.frame with the results will be returned.
  }
  \item{details}{
If \code{TRUE} the design characteristics and the steps during
sample size calculations will be shown.\cr 
If not specified, the default value is \code{FALSE} for \code{prior.type != "both"}
and \code{TRUE} otherwise.
  }
}
\details{
The sample size is calculated based on iterative evaluation of expected power.
The starting value of the sample size search is taken from a large sample 
approximation if \code{prior.type = "CV"}. Otherwise, an empirical start value is
obtained. Note that in case of \code{prior.type = "both"} the calculation may 
still take several seconds.

Note also that the expected power is always bounded above by the so-called
probability of technical success (\acronym{PTS}) which may be a value less than 1. 
Therefore, it may be possible that it is either not possible to calculate the 
required sample size at all or that the sample size gets very large 
if the given targetpower is less but close to the \acronym{PTS}.\cr\cr
The estimated sample size gives always the \emph{total} number of subjects (not subject/sequence in crossovers or subjects/group in parallel designs -- like in some other software packages).
}
\value{
A data.frame with the input values and the result of the sample size estimation.\cr
The \code{Sample size} column contains the total sample size in case of all
designs implemented.
}
\references{
Grieve AP. \emph{Confidence Intervals and Sample Sizes.} Biometrics. 1991;47:1597--603. \doi{10.2307/2532411}

O\enc{’}{'}Hagan, Stevens, JW, Campell MJ. \emph{Assurance in Clinical Trial Design.} Pharm Stat. 2005;4:187--201. \doi{10.1002/pst.175}

Julious SA, Owen RJ. \emph{Sample size calculations for clinical studies allowing for 
uncertainty in variance.} Pharm Stat. 2006;5:29--37. \doi{10.1002/pst.197}

Julious SA. \emph{Sample sizes for Clinical Trials.} Boca Raton: CRC Press / Chapman & Hall; 2010.

Bertsche A, Nehmitz G, Beyersmann J, Grieve AP. \emph{The predictive distribution of the residual variability in
the linear-fixed effects model for clinical cross-over trials.} Biom J. 2016;58(4):797--809. \doi{10.1002/bimj.201500245}

Box GEP, Tiao GC. \emph{Bayesian Inference in Statistical Analysis.} Boston: Addison-Wesley; 1992.

Held L, Sabanes Bove D. \emph{Applied Statistical Inference. Likelihood and Bayes.} Berlin, Heidelberg: Springer; 2014. \doi{10.1007/978-3-642-37887-4}

Senn S. \emph{Cross-over Trials in Clinical Research.} Chichester: John Wiley & Sons; 2\ifelse{html}{\out{<sup>nd</sup>}}{\eqn{^{nd}}} edition 2002.

Zierhut ML, Bycott P, Gibbs MA, Smith BP, Vicini P. \emph{Ignorance is not bliss: Statistical power is not probability of trial success.} Clin Pharmacol Ther. 2015;99:356--9. \doi{10.1002/cpt.257}
}
\author{
B. Lang, D. Labes
}
\seealso{
 \code{\link{exppower.TOST}, \link{known.designs}, \link{sampleN.TOST}}
}
\examples{
# Classical 2x2 cross-over, target power = 80\%,
# BE limits 80 ... 125\%, assumed true BE ratio = 95\%,
# intra-subject CV=30\% estimated from prior 2x2 trial 
# with m = 30 subjects
expsampleN.TOST(CV=0.3, prior.parm = list(m = 30, design = "2x2"))
# -> gives n = 42 with achieved expected power 0.806262
# Compare this to the usual sample size with CV assumed known ('carved in stone')
sampleN.TOST(CV=0.3)
# -> gives n = 40 subjects
# Compare this to the case where uncertainty is accounted for CV and theta0
# Not run due to timing policy of CRAN - may run several seconds
\donttest{
expsampleN.TOST(CV=0.3, prior.parm = list(m = 30, design = "2x2"), 
                prior.type = "both")}
# -> gives n = 72 subjects

# More than one CV with corresponding degrees of freedom 
# other settings as above in first example
CVs <- c(0.25, 0.3)
dfs <- c(22, 10)
expsampleN.TOST(CV=CVs, prior.parm = list(df = dfs))
# -> gives a pooled CV=0.2664927 with df=32
# and a sample size n=34 with achieved expected power 0.812653 exact
# achieved expected power 0.815019 approximate acc. Julious
}
