# ===========================================
# File: R/dimv.R
# Core R-native DIMV implementation
# - Conditional ridge-based imputation
# - Deterministic + Multiple-imputation variant
# ===========================================

#' Train DIMV Imputer (Ridge-Conditional)
#'
#' Trains a conditional ridge regression model for each variable to impute
#' missing values. Missing entries are first initialized using column means,
#' and then iteratively refined until convergence. Residual variances are
#' stored for use in multiple imputation.
#'
#' @param X A numeric data frame or matrix (may contain NA).
#' @param lambda Ridge penalty (default = 1.0).
#' @param maxit Maximum number of iterations (default = 50).
#' @param tol Convergence threshold for RMSE of updated missing values (default = 1e-4).
#' @param verbose If TRUE, prints iteration logs.
#'
#' @return An object of class `dimv_imputer` containing:
#'   - fitted conditional models (per column)
#'   - column means used for initialization
#'   - residual variances
#'   - convergence metadata
#'
#' @examples
#' set.seed(1)
#' X <- data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50))
#' # Inject missing values safely (5 per column)
#' for (j in seq_along(X)) {
#'   X[sample(1:nrow(X), 5), j] <- NA
#' }
#' imp <- dimv_train(X)
#' print(imp)
#' @export
#' @importFrom MASS ginv
dimv_train <- function(X, lambda = 1.0, maxit = 50, tol = 1e-4, verbose = FALSE) {
  # Ensure all variables are numeric
  X <- as.data.frame(lapply(X, as.numeric))
  n <- nrow(X); p <- ncol(X)
  coln <- colnames(X)
  if (is.null(coln)) coln <- paste0("V", seq_len(p))
  colnames(X) <- coln
  
  # Initial mean imputation (warm start)
  col_means <- sapply(X, function(v) mean(v, na.rm = TRUE))
  X_imp <- as.matrix(X)
  nas <- is.na(X_imp)
  if (any(nas)) X_imp[nas] <- rep(col_means, each = n)[nas]
  
  # Internal helper: ridge regression for one variable conditional on the others
  ridge_fit <- function(Xmat, y, lambda) {
    Z <- cbind(1, Xmat)                      # add intercept
    D <- diag(c(0, rep(1, ncol(Xmat))))      # ridge penalty (no penalty for intercept)
    A <- crossprod(Z) + lambda * D
    b <- crossprod(Z, y)
    coef <- tryCatch(solve(A, b), error = function(e) MASS::ginv(A) %*% b)
    list(intercept = coef[1], beta = coef[-1])
  }
  
  prev <- X_imp
  models <- vector("list", length = p); names(models) <- coln
  it <- 0
  
  for (it in seq_len(maxit)) {
    for (j in seq_len(p)) {
      # Fit conditional model for column j
      yvec <- prev[, j]
      Xo <- prev[, -j, drop = FALSE]
      fit <- ridge_fit(Xo, yvec, lambda)
      models[[j]] <- fit
      
      # Update only the originally missing values in column j
      na_idx <- which(is.na(as.matrix(X)[, j]))
      if (length(na_idx)) {
        preds <- fit$intercept + as.numeric(Xo[na_idx, , drop = FALSE] %*% fit$beta)
        prev[na_idx, j] <- preds
      }
    }
    
    # Check convergence only on originally missing entries
    if (any(nas)) {
      diff <- sqrt(mean((prev[nas] - X_imp[nas])^2))
      if (is.na(diff) || is.nan(diff)) diff <- 0
      if (verbose) cat("iter", it, "diff =", diff, "\n")
      X_imp <- prev
      if (diff < tol) break
    } else {
      if (verbose) cat("iter", it, ": no missing; stopping.\n")
      break
    }
  }
  
  # Compute residual variance for each conditional model
  resid_var <- numeric(p)
  for (j in seq_len(p)) {
    yvec <- X_imp[, j]
    Xo <- X_imp[, -j, drop = FALSE]
    predj <- models[[j]]$intercept + Xo %*% models[[j]]$beta
    resid_var[j] <- mean((yvec - predj)^2)
  }
  
  structure(list(models   = models,
                 col_means = col_means,
                 resid_var = resid_var,
                 lambda    = lambda,
                 iters     = it,
                 colnames  = coln),
            class = "dimv_imputer")
}


#' Impute New Data with a Fitted DIMV Imputer
#'
#' Applies a trained DIMV imputer to new (or the same) data. Iteratively imputes
#' missing values until convergence, using previously learned conditional models.
#'
#' @param imputer A `dimv_imputer` object returned by `dimv_train()`.
#' @param X_new A numeric data frame or matrix with missing values.
#' @param maxit Maximum number of update iterations.
#' @param tol Convergence threshold.
#' @param verbose If TRUE, prints iteration logs.
#'
#' @return A fully imputed data frame.
#' @examples
#' # After running dimv_train():
#' # dimv_impute_new(imp, X)
#' @export
dimv_impute_new <- function(imputer, X_new, maxit = 50, tol = 1e-4, verbose = FALSE) {
  if (!inherits(imputer, "dimv_imputer")) stop("imputer must be result of dimv_train()")
  X_new <- as.data.frame(lapply(X_new, as.numeric))
  X_new <- X_new[, imputer$colnames, drop = FALSE]
  
  # Mean initialization for new data
  X_imp <- as.matrix(X_new)
  nas <- is.na(X_imp)
  for (j in seq_along(imputer$col_means)) X_imp[is.na(X_imp[, j]), j] <- imputer$col_means[j]
  prev <- X_imp
  
  for (it in seq_len(maxit)) {
    for (j in seq_along(imputer$models)) {
      mod <- imputer$models[[j]]
      na_idx <- which(nas[, j])
      if (length(na_idx)) {
        preds <- mod$intercept + as.numeric(prev[na_idx, -j, drop = FALSE] %*% mod$beta)
        prev[na_idx, j] <- preds
      }
    }
    
    if (any(nas)) {
      diff <- sqrt(mean((prev[nas] - X_imp[nas])^2))
      if (is.na(diff) || is.nan(diff)) diff <- 0
      if (verbose) cat("impute iter", it, "diff =", diff, "\n")
      X_imp <- prev
      if (diff < tol) break
    } else {
      if (verbose) cat("no missing to impute; stopping.\n")
      break
    }
  }
  
  as.data.frame(X_imp)
}


#' Multiple Imputation via DIMV (Rubin-ready)
#'
#' Generates \code{m} multiply-imputed datasets using the DIMV conditional
#' expectations as the mean structure, and adds Gaussian noise based on
#' the learned residual variances. Only originally missing entries receive
#' stochastic noise, ensuring valid Rubin-style variance pooling.
#'
#' @param imputer A `dimv_imputer` model fitted by `dimv_train()`.
#' @param X_new Data containing missing values.
#' @param m Number of imputations to generate (default = 5).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A list of \code{m} fully-imputed data frames.
#' @examples
#' # imp_list <- dimv_impute_multiple(imp, X, m = 5, seed = 123)
#' @export
#' @importFrom stats rnorm
dimv_impute_multiple <- function(imputer, X_new, m = 5, seed = NULL) {
  if (!inherits(imputer, "dimv_imputer")) stop("imputer must be result of dimv_train()")
  if (!is.null(seed)) set.seed(seed)
  
  X_new <- as.data.frame(lapply(X_new, as.numeric))
  coln <- imputer$colnames
  X_new <- X_new[, coln, drop = FALSE]
  
  # First generate the deterministic DIMV imputation
  base_imp <- dimv_impute_new(imputer, X_new)
  
  res_list <- vector("list", m)
  for (k in seq_len(m)) {
    Xk <- as.matrix(base_imp)
    nas <- is.na(as.matrix(X_new))  # noise only added to originally missing cells
    
    for (j in seq_len(ncol(X_new))) {
      if (any(nas[, j])) {
        sigma <- sqrt(max(0, imputer$resid_var[j]))
        Xk[nas[, j], j] <- Xk[nas[, j], j] + stats::rnorm(sum(nas[, j]), mean = 0, sd = sigma)
      }
    }
    
    res_list[[k]] <- as.data.frame(Xk)
  }
  
  res_list
}


#' @exportS3Method print dimv_imputer
print.dimv_imputer <- function(x, ...) {
  cat("R-native DIMV imputer (ridge-conditional)\n",
      "lambda =", x$lambda, "| iterations =", x$iters, "\n")
}
