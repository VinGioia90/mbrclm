#' mbrclm: get the summary
#'
#' @description  Such a function allows for obtaining the results of the fit
#' @param object object of class mbrclm
#' @param ... further arguments
#'
#' @return The summary is printed
#' @export summary.mbrlcm
#' @export
#'
#' @name summary
#' @examples
summary.mbrclm <- function(object,   ...)
{
  print.summary.mbrclm <- function(x, digits = max(3, getOption("digits") - 3), ...)
  {
    #cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

    if(!x$convergence) {
      cat("model did not converge\n")
    }


    cat(paste("\nCoefficients (with ", x$link, " link):\n", sep = ""))
    printCoefmat(x$coefficients, digits = digits, signif.legend = FALSE)

    if(getOption("show.signif.stars") & any( x$coefficients[, 4L] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

    cat("\nType of estimator:", x$typeestimator, switch(x$typeestimator,
                                                        "AS_ml" = "(maximum likelihood)",
                                                        "AS_mean" = "(mean bias-reduced)",
                                                        "AS_median" = "(median bias-reduced)"))
    cat(paste("\nNumber of iterations in the quasi-Fisher scoring:", x$iterations, "\n"))
    cat(paste("\nConvergence status:", x$convergence, "\n"))
    invisible(x)
  }
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
  print.summary.mbrclm(object)
  return(invisible(object))
}
