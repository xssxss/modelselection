#' Title: Fitting linear models for the best model
#'
#' Description: `lm.best` is used to fit linear model for the best model provided by `ModelSelect`.
#'
#' @param object the model selection result from `ModelSelect`.
#' @param method the criteria to do model select.
#'               `method = "model"` selects the best model by the highest posterior probabilities.
#'               `method = "variable"` selects the variables in the best model by the posterior inclusion probabilities which are larger than the threshold.
#' @param threshold The threshold for variable selection. The variables with posterior inclusion probability larger than the threshold are selected in the best model. The default is 0.95.
#' @param x,y logicals. If `TRUE` the corresponding components (the best model predictor matrix, the response) of the fit are returned.
#' @return An object of class `"lm"`, which is a list containing the following components:
#' \describe{
#'   \item{\code{coefficients}}{A named vector of coefficients.}
#'   \item{\code{residuals}}{The residuals, that is the response minus the fitted values.}
#'   \item{\code{fitted.values}}{The fitted mean values.}
#'   \item{\code{rank}}{The numeric rank of the fitted linear model.}
#'   \item{\code{df.residual}}{The residual degrees of freedom.}
#'   \item{\code{call}}{The matched call.}
#'   \item{\code{terms}}{The `terms` object used.}
#'   \item{\code{model}}{(If requested) the model frame used.}
#'   \item{\code{qr}}{(If requested) the QR decomposition of the design matrix.}
#'   \item{\code{xlevels}}{(If the model formula includes factors) a record of the levels of the factors.}
#'   \item{\code{contrasts}}{(If the model formula includes factors) the contrasts used.}
#'   \item{\code{offset}}{The offset used.}
#'   \item{\code{threshold}}{the threshold used for method = "variable".}
#' }

#' @export
lm.best <- function(object, method = "model", threshold = 0.95, x = FALSE, y = FALSE){

  if(is.null(object)){
    stop("please provide object.")
  }

  X <- object$predictor
  Y <- object$response
  p <- ncol(X)
  n <- length(Y)
  response_name <- attr(Y, "id")
  if(method == "model"){
    index <- which(object[[method]][1,1:p]==1)
    if(length(index) == 0){
      fit <- stats::lm(Y~1, x=x, y=y)
      fit$call <- call("lm", stats::as.formula(paste0(response_name, "~1")))
    }else{
      Xsub <- X[,index]
      fit <- stats::lm(Y~Xsub, x=x, y=y)
      fit$call <- call("lm", stats::as.formula(paste0(response_name,"~", paste(colnames(Xsub), collapse ="+"))))
      names(fit$coefficients)[2:(length(index)+1)] <- colnames(Xsub)
    }
  }

  if(method == "variable"){
    index <- which(object[[method]][,"PIP"]>threshold)
    if(length(index) == 0){
      fit <- stats::lm(Y~1, x=x, y=y)
    }else{
      Xsub <- X[,index]
      fit <- stats::lm(Y~Xsub, x=x, y=y)
      fit$call <- call("lm", stats::as.formula(paste0(response_name,"~", paste(colnames(Xsub), collapse ="+"))))
      names(fit$coefficients)[2:(length(index)+1)] <- colnames(Xsub)
    }
    fit[["threshold"]] <- threshold
  }

  cat("Call:\n")
  print(fit$call)

  cat("\nCoefficients:\n")
  print(fit$coefficients)
  invisible(fit)

}


