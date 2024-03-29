\encoding{utf-8}
\name{pa.NTID}
\alias{pa.NTID}
\alias{pa.NTIDFDA}
\title{
Power analysis for scaled ABE for NTIDs
}
\description{
An analysis tool for exploration/visualization of the impact of expected values
(CV, theta0, reduced sample size due to drop-outs) on power of BE decision via
scABE for narrow therapeutic drugs (NTIDs) if these values deviate from the ones
assumed in planning the sample size of the study.\cr
The only implemented design is the full replicate design \code{"2x2x4"} according to the
FDA\enc{’}{'}s warfarin guidance.
}
\usage{
pa.NTID(CV, theta0 = 0.975, targetpower = 0.8, minpower = 0.7, ...)
}
\arguments{
  \item{CV}{
Coefficient of variation of the intra-subject variabilities of Test and Reference
as ratio (not percent).\cr
Here only the case \code{CVwT == CVwR} is implemented, \emph{i.e.}, CV has to be a scalar (\code{length(CV) == 1}).
}
  \item{theta0}{
\sQuote{True} or assumed T/R ratio. Often named GMR.\cr
Must be given as ratio. Defaults here to 0.975.
}
  \item{targetpower}{
Power to achieve at least in sample size estimation. Must be >0 and <1.\cr
Typical values are 0.8 or 0.9. Defaults to 0.8.\cr
Note that targetpower < 0.5 doesn\enc{’}{'}t make much sense.
}
  \item{minpower}{
Minimum acceptable power to have if deviating from assumptions for sample size plan.\cr
Has to be lower than \code{targetpower}. Defaults to 0.7.\cr
\code{minpower} < 0.5 doesn\enc{’}{'}t make much sense.
}
  \item{\dots}{
More arguments to pass to \code{power.NTID()}.\cr
F.i., \code{alpha}, \code{theta1}, \code{theta2} or \code{nsims} if other values
than the defaults for these arguments are needed. \cr
See man page of \code{power.NTID()}.
}
}
\details{
Power calculations are done via \code{power.NTID()} and
calculations of \code{CV} and \code{theta0} which result in \code{minpower} are obtained via \code{uniroot()}.\cr
While one of the parameters (\code{CV}, \code{theta0}, \code{n}) is varied, the respective two others are
kept constant. The tool shows the relative impact of single parameters on power.\cr
The tool takes a minimum of 12 subjects into account as demanded in most BE guidances.
However, it should be kept in mind that the FDA requires at least 24 subjects to be enrolled
in studies intended for reference-scaling.\cr\cr
It should be kept in mind that this is \bold{not} a substitute for the \sQuote{Sensitivity Analysis}
recommended in ICH-E9. In a real study a combination of all effects occurs simultaneously.
It is up to \emph{you} to decide on reasonable combinations and analyze their respective power.
}
\value{
Returns a list with class \code{'pwrA'} with the components
  \item{plan}{A data.frame with the result of the sample size estimation.
              See output of \code{sampleN.NTID()}.}
  \item{paCV}{A data.frame with value pairs CV, pwr for impact of
              deviations from CV.}
  \item{paGMR}{A data.frame with value pairs theta0, pwr for impact of
               deviations from theta0 (GMR).}
  \item{paN}{A data.frame with value pairs N, pwr for impact of
             deviations from planned N (dropouts).}
  \item{method}{Method of BE decision. Here "NTID".}
  \item{regulator}{Here "FDA".}
  \item{minpower}{Minimum acceptable power from the call of the function.}

The class \code{'pwrA'} has the S3 methods \code{print()} and \code{plot()}.
See \code{\link{pa.ABE}} for usage.
}
\references{

Food and Drug Administration, Office of Generic Drugs (OGD). \emph{Draft Guidance on Warfarin Sodium.} Recommended Dec 2012. \href{https://www.accessdata.fda.gov/drugsatfda_docs/psg/Warfarin_Sodium_tab_09218_RC12-12.pdf}{download}

Food and Drug Administration, Center for Drug Evaluation and Research (CDER). \emph{Draft Guidance for Industry. Bioequivalence Studies with Pharmacokinetic Endpoints for Drugs Submitted Under an ANDA.} August 2021. \href{https://www.fda.gov/media/87219/download}{download}

Yu LX, Jiang W, Zhang X, Lionberger R, Makhlouf F, Schuirmann DJ, Muldowney L, Chen ML, Davit B, Conner D, Woodcock J. \emph{Novel bioequivalence approach for narrow therapeutic index drugs.} Clin Pharmacol Ther. 2015;97(3):286--91. \doi{10.1002/cpt.28}

Jiang W, Makhlouf F, Schuirmann DJ, Zhang X, Zheng N, Conner D, Yu LX, Lionberger R. \emph{A Bioequivalence Approach for Generic Narrow Therapeutic Index Drugs: Evaluation of the Reference-Scaled Approach and Variability Comparison Criterion.} AAPS J. 2015;17(4):891--901. \doi{10.1208/s12248-015-9753-5}

\enc{Endrényi}{Endrenyi} L, \enc{Tóthfalusi}{Tothfalusi} L. \emph{Determination of Bioequivalence for Drugs with Narrow Therapeutic Index: Reduction of the Regulatory Burden.} J Pharm Pharm Sci. 2013;16(5):676--82. \href{https://journals.library.ualberta.ca/jpps/index.php/JPPS/article/download/20900/15927/0}{open access}
}
\author{
D. Labes according to code by H. \enc{Schütz}{Schuetz} for \code{pa.ABE()} and \code{pa.scABE()}.
}
\note{
The code for impact of deviations from planned sample size tries
to keep the degree of imbalance as low as possible between (sequence) groups.
This results in a lesser decrease of power than more extreme dropout-patterns.
}
\section{Warning }{
Be extremly carefull if your sample size plan has extremly small CV near or
below 0.05 (5\%). Adapt in that case your expected true ratio (\code{theta0})
to values nearer to 1 to not run into errors and/or long execution times.
}
\seealso{
  \code{\link{power.NTID}, \link{print.pwrA}, \link{plot.pwrA}, \link{pa.ABE}, \link{pa.scABE}}
}
\examples{
# using the defaults:
# targetpower=0.8, minpower=0.7, theta0/GMR=0.975
# BE margins from defaults of sampleN.NTID() 0.9002 ... 1.1108
# 1E5 sims in power.NTID()
# not run due to timing policy of CRAN for examples
# may run some ten seconds or more
\donttest{
plot(pa.NTID(CV=0.1))}
}
