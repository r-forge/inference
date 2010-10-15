\name{inference-class}
\alias{inference-class}
\alias{inference}
\title{An S4 class that stores inferential values of a fitted model object.}
\description{An S4 class that stores inferential values of a fitted model object.}
\details{An S4 class that inherits from the \code{matrix} class.  Rows
correspond to different coefficients and columns consist of point
estimates (point.est), confidence intervals (ci.lo and ci.hi),
p-values (p.value), and sample size (n).}
\docType{class}
\section{Slots}{\describe{\item{\code{.Data}:}{}\item{\code{model}:}{(\code{\link{character}}) Class of model fit; object of class \code{character},
such as "lm".}\item{\code{sample.size}:}{(\code{\link{numeric}}) Sample size used in model fit; object of class
\code{numeric}.}\item{\code{robust.se}:}{(\code{\link{logical}}) Boolean indicator whether robust standard errors were
used; object of class \code{logical}.}\item{\code{two.sided}:}{(\code{\link{logical}}) Boolean indicator whether p-values corresond to a
two-sided test or one-sided; object of class \code{logical}.}\item{\code{ci.level}:}{(\code{\link{numeric}}) Confidence level; object of class \code{numeric}.}\item{\code{scale}:}{(\code{\link{character}}) Scale of point estimates; object of class
\code{character}; defaults to "beta".}\item{\code{others}:}{(\code{\link{list}}) List containing other information about the model;
eg, summary of cluster size for \code{gee} and \code{lme} objects;
number of events for \code{coxph} objects.}}}
\section{Extends}{\code{\linkS4class{matrix}}}

