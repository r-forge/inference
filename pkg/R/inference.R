##' Extract inferential information from different statistical models.
##'
##' This package provides functions to extract point estimates, standard
##' errors, confidence intervals, p-values, and sample size from a fitted
##' model in a matrix-like object.  The purpose is to have all inferential
##' numbers be readily accessible, especially in the construction
##' of summary tables (R -> LaTeX -> html -> Word) for publication and
##' collaboration.
##'
##' @name inference-package
##' @rdname inference-package
##' @aliases inference
##' @docType package
##' @author Vinh Nguyen \email{vinhdizzo at gmail dot com}
##' @keywords statistical models point estimates confidence intervals
##' p-values sample size inference
##' @import methods
##' @examples
##' infer(lm(rnorm(100) ~ runif(100)))
NULL

##' An S4 class that stores inferential values of a fitted model object.
##'
##' An S4 class that inherits from the \code{matrix} class.  Rows
##' correspond to different coefficients and columns consist of point
##' estimates (point.est), confidence intervals (ci.lo and ci.hi),
##' p-values (p.value), and sample size (n).
##'
##' @rdname inference-class
##' @docType class
##' @slot .Data Object of class \code{matrix}.
##' @slot model Class of model fit; object of class \code{character},
##' such as "lm".
##' @slot sample.size Sample size used in model fit; object of class
##' \code{numeric}.
##' @slot robust.se Boolean indicator whether robust standard errors were
##' used; object of class \code{logical}.
##' @slot two.sided Boolean indicator whether p-values corresond to a
##' two-sided test or one-sided; object of class \code{logical}.
##' @slot ci.level Confidence level; object of class \code{numeric}.
##' @slot others List containing other information about the model;
##' eg, summary of cluster size for \code{gee} and \code{lme} objects;
##' number of events for \code{coxph} objects.
##' @exportClass inference
setClass(Class="inference"
         , representation=representation(model="character"
             , sample.size="numeric"
             , robust.se="logical"
             , two.sided="logical"
             , ci.level="numeric"
             , others="list")
         , contains=c("matrix"))

##' Show/print \code{inference} object.
##'
##' \code{show} method for objects made using the \code{infer} function.
##' @rdname show,inference-method
##' @aliases print.inference
##' @param object \code{inference} object.
##' @return Nothing.
##' @author Vinh Nguyen
setMethod("show", "inference", function(object){ print(slot(object, ".Data"))})

##' Inference for fitted model objects.
##'
##' Extract point estimates, standard errors, confidence intervals,
##' p-values, and sample size.
##' @rdname infer,-methods
##' @aliases infer infer,-method infer,lm-method infer,glm-method
##' @docType methods
##' @usage infer(fitobj, vars=NULL, robust.se=TRUE, two.sided=TRUE
##' , ci.level=0.95, ...)
##' @param fitobj Fitted model object, such as those of class \code{lm}.
##' @param vars Vector of variable names to obtain inference information
##' for.  Defaults to \code{NULL} which corresponds to all variables
##' in the fitted model.
##' @param robust.se Boolean indicator for whether robust standard
##' errors should be use.  Defaults to \code{TRUE}.
##' @param two.sided Boolean indicator for whether p-values should
##' correspond to a two-sided test or one-sided.  Defaults to
##' \code{TRUE}.
##' @param ci.level Confidence level.  Defaults to 0.95.
##' @param ... Not used.
##' @return S4 \code{inference} object.
##' @examples
##' infer(lm(rnorm(100) ~ runif(100)))
##' @exportMethod infer
##' @author Vinh Nguyen
setGeneric(name="infer", function(fitobj, ...) standardGeneric("infer"))

##' @nord
setMethod("infer", signature(fitobj="lm"), function(fitobj, vars=NULL, robust.se=TRUE, two.sided=TRUE, ci.level=0.95)
{
  if(is.null(vars)) vars <- names(coef(fitobj))
  point.est <- coef(fitobj)[vars]
  if(robust.se){
    require(sandwich)
    se <- sqrt(diag(sandwich(fitobj))[vars])
  } else{
    se <- sqrt(diag(vcov(fitobj))[vars])
  }
  p.value <- 1 - pnorm(abs(point.est/se))
  if(two.sided) p.value <- 2*p.value
  ci.lo <- point.est - abs(qnorm((1-ci.level)/2)) * se
  ci.hi <- point.est + abs(qnorm((1-ci.level)/2)) * se
  n <- length(fitobj$residuals)
  ##return(cbind(point.est, se, p.value, ci.lo, ci.hi, n))
  rslt <- new("inference", cbind(point.est, se, p.value, ci.lo, ci.hi, n), model=class(fitobj), sample.size=n, robust.se=robust.se, two.sided=two.sided, ci.level=ci.level, others=list("empty"))
  return(rslt)
})

##' @nord
setMethod("infer", signature(fitobj="glm"), function(fitobj, vars=NULL, robust.se=TRUE, two.sided=TRUE, ci.level=0.95, transform=NULL)
{
  if(is.null(vars)) vars <- names(coef(fitobj))
  point.est <- coef(fitobj)[vars]
  if(robust.se){
    require(sandwich)
    se <- sqrt(diag(sandwich(fitobj))[vars])
  } else{
    se <- sqrt(diag(vcov(fitobj))[vars])
  }
  p.value <- 1 - pnorm(abs(point.est/se))
  if(two.sided) p.value <- 2*p.value
  ci.lo <- point.est - abs(qnorm((1-ci.level)/2)) * se
  ci.hi <- point.est + abs(qnorm((1-ci.level)/2)) * se
  n <- length(fitobj$residuals)
  ##return(cbind(point.est, se, p.value, ci.lo, ci.hi, n))
  rslt <- new("inference", cbind(point.est, se, p.value, ci.lo, ci.hi, n), model=class(fitobj), sample.size=n, robust.se=robust.se, two.sided=two.sided, ci.level=ci.level, others=list("empty"))
  return(rslt)
})

##' @nord
setMethod("infer", signature(fitobj="gee"), function(fitobj, vars=NULL, robust.se=TRUE, two.sided=TRUE, ci.level=0.95, transform=NULL)
{
  if(is.null(vars)) vars <- names(coef(fitobj))
  point.est <- coef(fitobj)[vars]
  if(robust.se){
    ##require(sandwich)
    se <- sqrt(diag(fitobj$robust.variance)[vars])
  } else{
    se <- sqrt(diag(fitobj$naive.variance)[vars])
  }
  p.value <- 1 - pnorm(abs(point.est/se))
  if(two.sided) p.value <- 2*p.value
  ci.lo <- point.est - abs(qnorm((1-ci.level)/2)) * se
  ci.hi <- point.est + abs(qnorm((1-ci.level)/2)) * se
  n <- length(unique(fitobj$id))
  nObs <- fitobj$nobs
  summaryClusters <- summary(tapply(fitobj$id, fitobj$id, length))
  ##return(cbind(point.est, se, p.value, ci.lo, ci.hi, n))
  rslt <- new("inference", cbind(point.est, se, p.value, ci.lo, ci.hi, n), model=class(fitobj), sample.size=n, robust.se=robust.se, two.sided=two.sided, ci.level=ci.level, others=list(nObs=nObs, summaryClusters=summaryClusters))
  return(rslt)
})

##' @nord
setMethod("infer", signature(fitobj="lme"), function(fitobj, vars=NULL, robust.se=FALSE, two.sided=TRUE, ci.level=0.95, transform=NULL)
{
  if(is.null(vars)) vars <- names(coef(fitobj))
  point.est <- fixed.effects(fitobj)[vars] ##coef(fitobj)[vars]
  if(robust.se){
    ##require(sandwich)
    ##se <- sqrt(diag(sandwich(fitobj))[vars])
    stop("Robust standard errors are not available with Linear Mixed Effects Models.")
  } else{
    se <- sqrt(diag(vcov(fitobj))[vars])
  }
  p.value <- 1 - pnorm(abs(point.est/se))
  if(two.sided) p.value <- 2*p.value
  ci.lo <- point.est - abs(qnorm((1-ci.level)/2)) * se
  ci.hi <- point.est + abs(qnorm((1-ci.level)/2)) * se
  n <- fitobj$dims$ngrps["Subject"]
  nObs <- fitobj$dims$N
  summaryClusters <- summary(tapply(fitobj$groups[, 1], fitobj$groups[, 1], length))
  ##return(cbind(point.est, se, p.value, ci.lo, ci.hi, n))
  rslt <- new("inference", cbind(point.est, se, p.value, ci.lo, ci.hi, n), model=class(fitobj), sample.size=n, robust.se=robust.se, two.sided=two.sided, ci.level=ci.level, others=list(nObs=nObs, summaryClusters=summaryClusters))
  return(rslt)
})

##' @nord
setMethod("infer", signature(fitobj="coxph"), function(fitobj, vars=NULL, robust.se=TRUE, two.sided=TRUE, ci.level=0.95, transform=NULL)
{
  if(is.null(vars)) vars <- names(coef(fitobj))
  point.est <- coef(fitobj)[vars]
  if(robust.se){
    require(sandwich)
    if(any(names(fitobj) == "naive.var")){
      se <- sqrt(diag(vcov(fitobj))[vars])
    } else{
      se <- sqrt(diag(sandwich(fitobj))[vars])
    }
  } else{
    if(any(names(fitobj) == "naive.var")){
      warning("Only robust standard errors available due to robust=TRUE specification in coxph.")
    }
    se <- sqrt(diag(vcov(fitobj))[vars])
  }
  p.value <- 1 - pnorm(abs(point.est/se))
  if(two.sided) p.value <- 2*p.value
  ci.lo <- point.est - abs(qnorm((1-ci.level)/2)) * se
  ci.hi <- point.est + abs(qnorm((1-ci.level)/2)) * se
  n <- length(fitobj$residuals)
  n.events <- sum(fitobj$y[, 2])
  ##return(cbind(point.est, se, p.value, ci.lo, ci.hi, n))
  rslt <- new("inference", cbind(point.est, se, p.value, ci.lo, ci.hi, n), model=class(fitobj), sample.size=n, robust.se=robust.se, two.sided=two.sided, ci.level=ci.level, others=list(n.events=n.events))
  return(rslt)
})

## library(roxygen)
## package.skeleton("inference2", code_files="inference.R", force=TRUE)
## roxygenize(package.dir="inference", roxygen.dir="inference", copy.package=FALSE, unlink.target=FALSE)
## roxygenize(package.dir="inference", roxygen.dir="inference", copy.package=FALSE, unlink.target=FALSE, use.Rd2=TRUE) ## will try to document all function declarations, even if I did not document it using roxygen
## roxygenize(package.dir="pkg", roxygen.dir="pkg", copy.package=FALSE, unlink.target=FALSE, use.Rd2=TRUE) ## use this for S4
## R CMD roxygen -d -s pkg
