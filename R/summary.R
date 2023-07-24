#' mbrclm: get the summary
#'
#' @description  Such a function allows for obtaining the results of the fit
#' @param object object of class mbrclm
#' @param ... further arguments
#'
#' @return Summary
#'
summary.mbrclm <- function(object,   ...)
{

  ## coefficient table

  cf <- as.vector(object$coefficients)
  se <- object$std.err
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object$coefficients <- cf

  ## number of iterations

  object$iterations <- object$iter

  ## type estimator
  object$typeestimator <- object$type


  ## return
  class(object) <- "summary.mbrclm"
  object
}
