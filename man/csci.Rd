\name{CSCI}
\alias{CSCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pointwise Confidence Intervals for Current Status Data
}
\description{
Calculates several different methods for getting pointwise confidence intervals for current st
}
\usage{
CSCI(C, D, times=NULL, type = c("VALID", "ABA", "LIKELIHOOD"), 
   conf.level = 0.95, control=controlCSCI())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{C}{
a vector of assessement times
}
  \item{D}{
a vector of indicators of event at or before the assessment time
}
  \item{times}{
a vector of times, t, to give the confidence interval for the event time distribution, F(t). If NULL then set to \code{sort(unique(C))}.
}
  \item{type}{
type of confidence interval, either "VALID", "ABA", or "LIKELIHOOD" (see details)
}
  \item{conf.level}{
confidence level for intervals (for \code{type="LIKELIHOOD"} only specific values are allowed, see note)
}
\item{control}{
list with parameters for algorithms, see \code{\link{controlCSCI}}
}
}
\details{
The function does three types of pointwise confidence intervals for the cumulative distribution function for the event time at the times specified by \code{times}. When \code{type="VALID"} the function gives a method that guarantees that the coverage will be at least nominal, but the confidence intervals are not ensured to be monotonic over the times of interest. When \code{type="ABA"} the function 
gives an approximate method that 
does not guarantee coverage, but has been shown by simulation to have good coverage for 
smoothly changing distributions,
and it does ensure monotonicity (see Kim, et al, 2020). 
When \code{type="LIKELIHOOD"} the function gives an asymptotic likelihood ratio test-based confidence interval that does not guarantee coverage
(Banerjee and Wellner, 2001).

}
\value{
A list with 2 objects:
\item{ciTable_all}{data.frame with NPMLE and associated confidence intervals for all possible time values (not output for type='LIKELIHOOD')}
\item{ciTable_times}{data.frame with NPMLE and assoicated confidence intervals for the values of 'times' argument}
%% ...
}
\references{
Banerjee, M. and J. A. Wellner (2001). Likelihood ratio tests for monotone functions.
Ann. Statist. 29 (6), 1699-1731.


Kim, S, Fay, MP, Proschan, MA (2020). Valid and Approximately Valid Confidence Intervals for Current Status Data.
(see \url{https://arxiv.org/abs/1805.06488}).
}
\author{
Sungwook Kim
}
\note{Because the likelihood ratio test goes to a non-standard asymptotic distribution, we do not calculate quantiles from that distribution, but take them from Table 2 of Banerjee and Wellner (2001). Because of this, when \code{type="LIKELIHOOD"} then 
conf.level must be one of 0.25,0.50,0.75,0.80,0.85,0.90,0.95, or 0.99.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
\examples{
data(hepABulg)
CSCI(C=hepABulg$age,D=hepABulg$testPos,type="VALID")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{survival}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{current status data}
% \concept{ ~cpt2 }
% Use only one concept per line.
