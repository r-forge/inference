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
##' @aliases inference
##' @docType package
##' @title Extract inferential information from different statistical
##' models.
##' @author Vinh Nguyen \email{vqnguyen@@uci.edu}
##' @keywords statistical models, point estimates, confidence intervals,
##' p-values, sample size, inference
##' @examples
##' infer(lm(rnorm(100) ~ runif(100)))
NULL

##' <description>
##'
##' <details>
##' @title bye
##' @param fitobj 
##' @return 1
##' @author Vinh Nguyen
setGeneric("genF", function(fitobj, ...) standardGeneric("genF"))
##' <description>
##'
##' <details>
##' @title bye
##' @param fitobj 
##' @return 1
##' @author Vinh Nguyen
setMethod("genF", signature(fitobj="lm"), function(fitobj){ print("lm") })
##' <description>
##'
##' <details>
##' @title hello
##' @param fitobj 
##' @return 1
##' @author Vinh Nguyen
setMethod("genF", signature(fitobj="glm"), function(fitobj){ print("glm") })



##' An S4 class that stores inferential values of a fitted model object.
##'
##' An S4 class that inherits from the \code{matrix} class.  Rows
##' correspond to different coefficients and columns consist of point
##' estimates (point.est), confidence intervals (ci.lo and ci.hi),
##' p-values (p.value), and sample size (n).
##'
##' @name inference-class
##' @slot .Data matrix consisting of inferential values
##' @slot sample.size Sample size used in model fit.
##' @slot robust.se Boolean indicator whether robust standard errors were
##' used.
##' @slot two.sided Boolean indicator whether p-values corresond to a
##' two-sided test or one-sided.
##' @slot ci.level Confidence level.
setClass(Class="inference"
         , representation=representation(model="character"
             , sample.size="numeric"
             , robust.se="logical"
             , two.sided="logical"
             , ci.level="numeric")
         , contains=c("matrix"))

##' Show/print inference object.
##'
##' show method for objects made using \code{infer()}
##' @title \code{show} method for \code{inference} class.
##' @param object 
##' @return something
##' @author Vinh Nguyen
setMethod("show", "inference", function(object){ print(slot(object, ".Data"))})

##' infer generic function.
##'
##' infer() used on fitted model objects, such as \code{lm} objects.
##' @title infer
##' @param fitobj 
##' @param ... 
##' @return something
##' @author Vinh Nguyen
setGeneric("infer", function(fitobj, ...) standardGeneric("infer"))

##' Inference for \code{lm} objects
##'
##' Extract point estimates, standard errors, confidence intervals,
##' p-values, and sample size
##' @title infer-lm
##' @name infer-methods
##' @rdname infer-methods
##' @param fitobj Fitted object.
##' @param vars Vector of variable names we wish to extract information
##' for.  Defaults to NULL which corresponds to all variables in the
##' fitted model.
##' @param robust.se Boolean indicator for whether robust standard
##' errors should be use.  Defaults to TRUE.
##' @param two.sided Boolean indicator for whether p-values should
##' correspond to a two-sided test or one-sided.  Defaults to TRUE.
##' @param ci.level Confidence level.  Defaults to 0.95.
##' @return S4 \code{inference} object.
##' @author Vinh Nguyen
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
  rslt <- new("inference", cbind(point.est, se, p.value, ci.lo, ci.hi, n), model=class(fitobj), sample.size=n, robust.se=robust.se, two.sided=two.sided, ci.level=ci.level)
  return(rslt)
})
##' <description>
##'
##' <details>
##' @title infer-glm
##' @name infer-methods
##' @rdname infer-methods
##' @param fitobj 
##' @param vars 
##' @param robust.se 
##' @param two.sided 
##' @param ci.level 
##' @param transform 
##' @return S4 \code{inference} object.
##' @author Vinh Nguyen
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
  rslt <- new("inference", cbind(point.est, se, p.value, ci.lo, ci.hi, n), model=class(fitobj), sample.size=n, robust.se=robust.se, two.sided=two.sided, ci.level=ci.level)
  return(rslt)
})

## library(roxygen)
## package.skeleton("inference2", code_files="inference.R", force=TRUE)
## roxygenize(package.dir="inference", roxygen.dir="inference", copy.package=FALSE, unlink.target=FALSE)
## roxygenize(package.dir="inference", roxygen.dir="inference", copy.package=FALSE, unlink.target=FALSE, use.Rd2=TRUE)
