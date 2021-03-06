\name{ex_ICLT}
\alias{ex_ICLT}
\docType{data}
\title{Virtual data set for interval censored and left truncated data}
\description{
The data \code{ex_ICLT} is a virtual data set created to show how to utilize the package. \code{ex_ICLT} is interval censored and left truncated data. See \code{ex_IC} for interval censored data.
}
\usage{data(ex_ICLT)}
\format{
  The data have the following columns:
  \describe{
    \item{\code{lowerIC}}{The lower limit of the censoring interval.}
    \item{\code{upperIC}}{The upper limit of the censoring interval.}
    \item{\code{trunc}}{The vector of left truncated points.}
    \item{\code{X1 - X6}}{The covariate vectors used for variable selection.}
  }
}
\examples{
library(ALassoSurvIC)
data(ex_IC) # 100 subjects and 6 covariates
print(ex_IC)
}

\keyword{datasets}
\seealso{\code{\link{ex_IC}}}

