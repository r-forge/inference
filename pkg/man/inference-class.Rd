\name{inference-class}
\alias{inference-class}
\alias{inference}
\title{An S4 class that stores inferential values of a fitted model object.}
\description{An S4 class that stores inferential values of a fitted model object.}
\details{An S4 class that inherits from the \code{\link[=structure-class]{matrix}}
class in the \code{methods} package; see \code{class?matrix}.  Rows
correspond to different coefficients and columns consist of point
estimates (point.est), confidence intervals (ci.lo and ci.hi),
p-values (p.value), and sample size (n).}
\docType{class}
\section{Slots}{\describe{\item{\code{.Data}:}{}\item{\code{model}:}{(\code{\link{character}}) String specifying class of model fit, such as "lm".}\item{\code{sample.size}:}{(\code{\link{numeric}}) Sample size used in model fit.}\item{\code{robust.se}:}{(\code{\link{logical}}) Boolean indicator whether robust standard errors were
used.}\item{\code{two.sided}:}{(\code{\link{logical}}) Boolean indicator whether p-values corresond to a
two-sided test or one-sided.}\item{\code{ci.level}:}{(\code{\link{numeric}}) Confidence level, e.g., 0.95.}\item{\code{scale}:}{(\code{\link{character}}) Scale of point estimates; defaults to "beta".}\item{\code{others}:}{(\code{\link{list}}) List containing other information about the model;
eg, summary of cluster size for \code{\link[gee]{gee}} and \code{\link[nlme]{lme}} objects;
number of events for \code{\link[survival]{coxph}} objects.}}}
\section{Extends}{\code{\linkS4class{matrix}}}

