\name{infer,-methods}
\alias{infer-methods}
\alias{infer}
\title{Inference for fitted model objects.}
\usage{infer(fitobj, vars=NULL, robust.se=TRUE, two.sided=TRUE
, ci.level=0.95, ...)}
\description{Inference for fitted model objects.}
\details{Extract point estimates, standard errors, confidence intervals,
p-values, and sample size.}
\alias{infer}
\alias{infer,-method}
\alias{infer,lm-method}
\alias{infer,glm-method}
\docType{methods}
\value{S4 \code{\linkS4class{inference}} object.}
\author{Vinh Nguyen}
\arguments{\item{fitobj}{Fitted model object, such as those of class \code{\link{lm}}.}
\item{vars}{Vector of variable names to obtain inference information
for.  Defaults to \code{NULL} which corresponds to all variables
in the fitted model.}
\item{robust.se}{Boolean indicator for whether robust standard
errors should be use.  Defaults to \code{TRUE}.}
\item{two.sided}{Boolean indicator for whether p-values should
correspond to a two-sided test or one-sided.  Defaults to
\code{TRUE}.}
\item{ci.level}{Confidence level.  Defaults to 0.95.}
\item{...}{Not used.}
}
\examples{infer(lm(rnorm(100) ~ runif(100)))}

