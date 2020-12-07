\name{controlCSCI}
\alias{controlCSCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for control parameters for algorithms used in CSCI.
}
\description{
Allows chainging of default parameters.
}
\usage{
controlCSCI(power = 2/3, 
   quan_p = c(0.25, 0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99), 
   xp_hat = c(0.06402, 0.28506, 0.80694, 0.98729, 1.22756, 1.60246, 2.26916, 3.8363), 
   intF = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{power}{
for defining m in the algorithm when \code{type='VALID'}: m=ceiling(n^power), where n=length(C)
}
  \item{quan_p}{
quantile associated with xp_hat, used when \code{type='LIKELIHOOD'}
}
  \item{xp_hat}{
estimated quantile of the distribution of the log likelihood ratio (see e.g., Table 2 of Banerjee and WWellner, 2001),
 used when \code{type='LIKELIHOOD'}
}
  \item{intF}{
numer of intervals to partition the F space (F=c(1:(intF-1)/intF)), 
used when \code{type='LIKELIHOOD'}
}
}
\details{
For \code{power}, see Kim, et al 2020. For details on the other values, see the code for the \code{type='LIKELIHOOD'} algorithm and Banerjee and Wellner, 2001. 
}
\value{
A list of the argument values. 
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Banerjee, M. and J. A. Wellner (2001). Likelihood ratio tests for monotone functions.
Ann. Statist. 29 (6), 1699-1731.


Kim, S, Fay, MP, Proschan, MA (2020). Valid and Approximately Valid Confidence Intervals for Current Status Data.
(see \url{https://arxiv.org/abs/1805.06488}).
}
%\author{%%  ~~who you are~~}
%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{%% ~~objects to See Also as \code{\link{help}}, ~~~}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
