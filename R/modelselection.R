#' Title: Model selection for linear models
#'
#' Description: use BIC to do model selection.
#'
#' @param formula an object of class "formula": a symbolic description of the model to be fitted.
#'                A typical model has the form `response ~ terms` where response is the (numeric) `response` vector and terms is a series of terms which specifies a linear predictor for `response`.
#'                A terms specification of the form `first + second` indicates all the terms in `first` together with all the terms in `second` with duplicates removed.
#'                A specification of the form `first:second` indicates the set of terms obtained by taking the interactions of all terms in `first` with all terms in `second`.
#'                The specification `first*second` indicates the cross of `first` and `second.` This is the same as `first + second + first:second`.
#' @param data an data frame containing the variables in the model.
#' @return `ModelSelect` returns a list containing the following components:
#'
#' \describe{
#'   \item{\code{model}}{A data frame of candidate models' BIC and posterior probabilities, sorted by decreasing posterior probability}
#'   \item{\code{variable}}{A data frame of candidate variables' posterior inclusion probabilities and decision}
#'   \item{\code{response}}{the response variable}
#'   \item{\code{predictor}}{All candidate predictor variables.}
#' }
#'
#' The function `lm.best` is used to obtain the linear fitting to the best model by posterior probability or by controlling variables' posterior inclusion probabilities.
#'
#' @export

ModelSelect <- function(formula, data){

  if(is.null(formula)){
    stop("provide formula.")
  }

  if(is.null(data)){
    stop("provide data.")
  }

  if(sum(is.na(data))){
    stop("missing data.")
  }

  terms_obj <- stats::terms(formula, data = data)
  terms_in_formula <- all.vars(terms_obj)
  response_name <- as.character(attr(terms_obj, "variables"))[attr(terms_obj, "response") + 1]

  missing_vars <- setdiff(terms_in_formula, names(data))

  if (length(missing_vars) > 0) {
    stop(paste("The following variables are missing from the data:", paste(missing_vars, collapse = ", ")))
  }

  model_frame <- stats::model.frame(formula = formula, data = data)
  Y <- stats::model.response(model_frame)
  X <- stats::model.matrix(object = formula, data = data)[,-1]


  if(!is.numeric(Y)){
    stop("Y should be numeric data.")
  }

  if(!is.numeric(X)){
    stop("X should be numeric data.")
  }


  #add formula example

  K <- ncol(X)
  if(K < 16){
    #low dimension
    #bic and posterior probability for each model
    model_bic_prob <- c()
    for (i in 0:(2^K-1)){
      Xc <- X[, which(binaryLogic::as.binary(i, n=K)), drop = FALSE]
      if(i == 0){
        fit <- stats::lm(Y~1)
        model_bic_prob <- rbind(model_bic_prob, c(binaryLogic::as.binary(i, n=K), stats::BIC(fit)))

      }else{

        sub_colname <- colnames(Xc)
        sub_var_name <- strsplit(sub_colname, split = ":")
        sub_var_len <- unlist(lapply(sub_var_name, length))

        if(!length(which(sub_var_len > 1))){
          fit <- stats::lm(Y~Xc)
          model_bic_prob <- rbind(model_bic_prob, c(binaryLogic::as.binary(i, n=K), stats::BIC(fit)))

        }else{
          sub_interaction <- sub_var_name[[which(sub_var_len > 1)]]
          sub_not_interaction <- sub_colname[which(sub_var_len == 1)]

          if(Reduce("&", sub_interaction %in% sub_not_interaction)){

            fit <- stats::lm(Y~Xc)
            model_bic_prob <- rbind(model_bic_prob, c(binaryLogic::as.binary(i, n=K), stats::BIC(fit)))

          }
        }
      }
    }

    bic <- model_bic_prob[,ncol(model_bic_prob)]
    bestmodel <- model_bic_prob[which.min(bic),-ncol(model_bic_prob)]
    postprob_temp <- exp(-0.5*(bic-min(bic)))
    postprob <- postprob_temp/sum(postprob_temp)
    model_bic_prob <- cbind(model_bic_prob, postprob)

    #postperior inclusion probability

    var_pip <- colnames(X)
    dat <- model_bic_prob[,1:K]
    pip <- colSums(dat * postprob)
    var_pip <- data.frame("Var"=var_pip, "bestmodel" = bestmodel,"PIP" = round(pip,3))

  }else{
    #high dimension
    # Do model search with genetic algorithm

    fitness_ftn <- function(string){
      if(sum(string) == 0){
        return(-stats::BIC(stats::lm(Y~1)))
      }else{

        model <- which(string==1)
        k <- length(model)
        Xsub <- X[,model,drop = FALSE]
        sub_colname <- colnames(Xsub)
        sub_var_name <- strsplit(sub_colname, split = ":")
        sub_var_len <- unlist(lapply(sub_var_name, length))
        sub_interaction <- sub_var_name[[which(sub_var_len > 1)]]
        sub_not_interaction <- sub_colname[which(sub_var_len == 1)]

        if(!length(which(sub_var_len > 1))){
          return(-stats::BIC(stats::lm(Y~Xsub)))
        }else{
          sub_interaction <- sub_var_name[[which(sub_var_len > 1)]]
          sub_not_interaction <- sub_colname[which(sub_var_len == 1)]
          if(sub_interaction %in% sub_not_interaction){

            return(-stats::BIC(stats::lm(Y~Xsub)))

          }else{
            return(-Inf)
          }
        }
      }
    }


    if(K > 99){
      suggestedsol <- diag(K)
      tmp_BIC <- vector()
      for(i in 1:K){
        model <- which(suggestedsol[i,]==1)
        k <- length(model)

        Xsub <- X[,model,drop = FALSE]
        tmp_BIC[i] <- -stats::BIC(stats::lm(Y~Xsub))
      }
      suggestedsol <- rbind(0,suggestedsol[order(tmp_BIC,decreasing = TRUE)[1:99],])
    }else{
      suggestedsol <- rbind(0,diag(K))
    }
    maxiterations = 2000
    runs_til_stop = 1000

    fitness_ftn <- memoise::memoise(fitness_ftn)
    ans <- GA::ga("binary", fitness = fitness_ftn,
                  nBits = K,maxiter = maxiterations,popSize = 100,
                  elitism = min(c(10,2^K)),run = runs_til_stop,suggestions = suggestedsol,monitor = TRUE)
    memoise::forget(fitness_ftn)
    dat <- ans@population
    dupes <- duplicated(dat)
    dat <- dat[!dupes,]
    ans@fitness <- ans@fitness[!dupes]
    bic <- -ans@fitness # BIC

    bestmodel <- dat[which.min(bic),] #ans@solution

    model_bic_prob <- cbind(dat,bic)
    postprob_temp <- exp(-0.5*(bic-min(bic)))
    postprob <- postprob_temp/sum(postprob_temp)
    model_bic_prob <- cbind(model_bic_prob, postprob)
    model_bic_prob <- model_bic_prob[-which(bic == Inf),]

    #postperior inclusion probability
    var_pip <- colnames(X)
    pip <- colSums(dat * postprob)
    var_pip <- data.frame("Var"=var_pip, "bestmodel" = bestmodel,"PIP" = round(pip,3))

  }
  model_bic_prob <- as.data.frame(model_bic_prob)
  if(is.null(colnames(X))){
    colnames(model_bic_prob) <- c(paste("X",1:K,sep=""),"BIC","posteriorprob")
  }else{
    colnames(model_bic_prob) <- c(colnames(X),"BIC","posteriorprob")

  }
  model_bic_prob <- model_bic_prob[order(model_bic_prob$posteriorprob, decreasing = T),]


  attr(model_bic_prob, "id") <- "model"
  attr(var_pip, "id") <- "variable"
  attr(Y, "id") <- response_name

  cat("The Best Model:\n", colnames(model_bic_prob)[1:K][model_bic_prob[1,1:K] == 1], "\n")
  cat("BIC:\n", model_bic_prob[1,K+1])
  invisible(list("model" = model_bic_prob, "variable" = var_pip, "response" = Y, "predictor" = X))

}



