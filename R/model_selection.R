
# Functions for performing model selection --------------------------------------------------------

#' @title Compute variance of ensemble mean
#' @param X n x m matrix of predictions for n data points across m models
#' @param S n x m matrix of estimate standard errors
#' @param w m-length vector of model weights
#'
#' @return A length n vector of ensemble mean standard errors
#' @export
weighted_se <- function(X, S, w) {
  
  var_x <- rep(0, nrow(X))
  
  for (i in 1:ncol(S)) {
    for (j in 1:ncol(S)) {
      var_ij <- w[i] * w[j] * cor(X[,i], X[,j], use = "complete.obs") * S[,i] * S[,j]
      var_x <- var_x + var_ij
    }
  }
  
  return(sqrt(var_x))
  
}

#' @title Obtain stacked weights
#' @param yi Length-n vector of observations
#' @param yhat n by m matrix of n predictions for m models
#'
#' @return A vector of length m corresponding to model weights
#' @export
solve_weights <- function(yi, yhat) {
  
  obj_fn <- function(ws, yi, yhat) {
    w <- plogis(ws)/sum(plogis(ws))
    W <- matrix(rep(w, nrow(yhat)), ncol = ncol(yhat), nrow = nrow(yhat), byrow = TRUE)
    sqrt(sum((yi - rowSums((W * yhat)))^2)/nrow(yhat))
  }
  
  opt <- optim(rep(0, ncol(yhat)), obj_fn, yi = yi, yhat = yhat)
  
  return(plogis(opt$par)/sum(plogis(opt$par)))
  
}

#' @title Simple function to linearly rescale variables to 0-1 interval
#' @param x Variable to rescale
#' @param na.rm Defaults to TRUE
#'
#' @return A numeric vector
#' @export
rescale_01 <- function(x, na.rm = TRUE) {
  (x - min(x, na.rm = na.rm)) / (max(x, na.rm = na.rm) - min(x, na.rm = na.rm))
}


#' @title complimentary log-log link / transform
#' @param x Numeric
#' @param inverse if TRUE, use inverse link
#'
#' @return A numeric vector
#' @export
cloglog <- function(x, inverse = FALSE) {
  if (!inverse) {
    1 - exp(-exp(x))
  } else {
    log(-log(1 - x))
  }
}


#' @title Wrapper function to extract AUC using dismo::evaluate
#' @param obs Observed presence-absence (0-1) data
#' @param fit Fitted probabilities
#'
#' @return AUC value
#' @export
AUC <- function(obs, fit) {
  p <- c(fit[obs == 1])
  a <- c(fit[obs == 0])
  dismo::evaluate(p, a)@auc
}