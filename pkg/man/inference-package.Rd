\name{inference-package}
\alias{inference-package}
\alias{inference-package}
\title{Extract inferential information from different statistical models.}
\description{Extract inferential information from different statistical models.}
\details{This package provides functions to extract point estimates, standard
errors, confidence intervals, p-values, and sample size from a fitted
model in a matrix-like object.  The purpose is to have all inferential
numbers be readily accessible, especially in the construction
of summary tables (R -> LaTeX -> html -> Word) for publication and
collaboration.}
\alias{inference}
\docType{package}
\author{Vinh Nguyen \email{vinhdizzo at gmail dot com}}
\keyword{statistical}
\keyword{models}
\keyword{point}
\keyword{estimates}
\keyword{confidence}
\keyword{intervals}
\keyword{p-values}
\keyword{sample}
\keyword{size}
\keyword{inference}
\examples{infer(lm(rnorm(100) ~ runif(100)))}

