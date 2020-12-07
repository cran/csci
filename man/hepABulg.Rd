\name{hepABulg}
\alias{hepABulg}
\docType{data}
\title{
Hepatitis A Data from Bulgaria
}
\description{
Hepatitis A data from Bulgaria, collected from school-children and blood donors by Prof. G. Frosner, Munich (from Keiding, 1991, Table 1).
}
\usage{data("hepABulg")}
\format{
  A data frame with 850 observations on the following 2 variables.
  \describe{
    \item{\code{age}}{a numeric vector}
    \item{\code{testPos}}{a numeric vector, Hepatitis A positive=1, or not=0}
  }
}
\details{
Each row in the data frame represents an individual and the age tested in years and the results of the hepatitis A test (1=positive, 0=negative). Ages of the individuals range from 1 to 86 years old.
}
\source{
Keiding, N (1991). Age-specific Incidence and Prevalence: a Statistical Perspective. JRSS A 154(3): 371-412 (Table 2).
}
%\references{%%  ~~ possibly secondary sources and usages ~~}
\examples{
data(hepABulg)
head(hepABulg)
}
\keyword{datasets}
