\name{infer,-methods}
\alias{infer-methods}
\alias{infer}
\title{infer}
\usage{infer(fitobj, vars, robust.se, two.sided=TRUE, ci.level=0.95, ...)}
\description{Inference for fitted model objects.}
\details{Extract point estimates, standard errors, confidence intervals,
p-values, and sample size.}
\alias{infer,lm-method}
\alias{infer,glm-method}
\alias{infer,coxph-method}
\alias{infer,gee-method}
\alias{infer,lme-method}
\alias{infer,mer-method}
\alias{infer,lm-method}
\alias{infer,lm,ANY,ANY,ANY,ANY-method}
\alias{infer,glm,ANY,ANY,ANY,ANY-method}
\alias{infer,coxph,ANY,ANY,ANY,ANY-method}
\alias{infer,gee,ANY,ANY,ANY,ANY-method}
\alias{infer,lme,ANY,ANY,ANY,ANY-method}
\alias{infer,mer,ANY,ANY,ANY,ANY-method}
\docType{methods}
\value{S4 \code{\linkS4class{inference}} object.}
\author{Vinh Nguyen}
\arguments{\item{fitobj}{Fitted model object, such as those of class \code{\link[stats]{lm}}.}
\item{vars}{Vector of variable names to obtain inference information
for.  If not specified, all variables in the fitted model will be used in the fitted model.}
\item{robust.se}{Boolean indicator for whether robust standard
errors should be use.  Defaults to \code{TRUE} for all models except \code{lme} and \code{mer}.}
\item{two.sided}{Boolean indicator for whether p-values should
correspond to a two-sided test or one-sided.  Defaults to
\code{TRUE}.}
\item{ci.level}{Confidence level.  Defaults to 0.95.}
\item{...}{Not used.}
}
\examples{infer(lm(rnorm(100) ~ runif(100)))}

