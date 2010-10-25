\name{transform.inference}
\alias{transform.inference}
\title{Transformation of point estimates}
\usage{transform.inference <- function(`_data`, f, f.prime, ...)}
\description{\code{\link[base]{transform}} method for class inference}
\details{Transform the point estimates, confidence intervals, and standard
errors based on the delta method.  This builds on the S3 generic
function \code{\link[base]{transform}} from the \code{base} package.

It can be used to get the hazard ratio scale in inference objects
created from \code{\link[survival]{coxph}} objects and the odds ratio scale from
logistic regression created from \code{\link[stats]{glm}}
(both using \code{f=exp, f.prime=exp}).}
\value{Object of class \code{\linkS4class{inference}}.}
\author{Vinh Nguyen}
\arguments{\item{`_data`}{Object of class \code{\linkS4class{inference}}.}
\item{f}{Function to transform the point estimates and confidence intervals; e.g., \code{\link[base]{exp}}.}
\item{f.prime}{Derivative of \code{f} in order to compute the standard
error of the transformed point estimates based on the delta method.}
\item{...}{Nothing.}
}

