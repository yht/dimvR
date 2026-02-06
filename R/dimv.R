# ===========================================
# File: R/dimv.R
# Core R-native DIMV implementation
# + Feature selection (optional; OFF by default, uses column INDEXES)
# + Adaptive regularization (optional; OFF by default)
# + Basic diagnostics helper (+ diff_history)
# Note: Feature selection & adaptive regularization are EXPERIMENTAL.
# ===========================================

#' @keywords internal
#' @noRd
# Simple feature selector based on absolute correlation (fallback)
feature_select_fn <- function(X, target_col, max_features = 5) {
  cor_vals <- suppressWarnings(stats::cor(X, use = "pairwise.complete.obs"))
  target_cor <- abs(cor_vals[, target_col])
  target_cor[target_col] <- NA
  best <- names(sort(target_cor, decreasing = TRUE))[1:max_features]
  best_idx <- match(best[!is.na(best)], colnames(X))
  best_idx[!is.na(best_idx)]
}

#' @keywords internal
#' @noRd
# Adaptive lambda based on condition number (ill-conditioning)
adaptive_lambda_fn <- function(X, base_lambda = 1.0) {
  k <- tryCatch(base::kappa(X), error = function(e) 1e6)
  base_lambda * log1p(k)
}

#' Train DIMV Imputer (Ridge-Conditional)
#'
#' Trains a conditional ridge regression model for each variable to impute
#' missing values. Missing entries are first initialized using column means,
#' and then iteratively refined until convergence. Optionally supports
#' per-variable feature selection (using column indexes) and adaptive
#' regularization (both disabled by default).
#'
#' @param X A numeric data frame or matrix (may contain \code{NA}).
#' @param lambda Base ridge penalty (default \code{1.0}).
#' @param maxit Maximum number of iterations (default \code{50}).
#' @param tol Convergence threshold for RMSE of updated missing values (default \code{1e-4}).
#' @param feature_select Logical; if \code{TRUE}, apply per-variable feature selection (default \code{FALSE}).
#' @param adaptive Logical; if \code{TRUE}, use adaptive lambda per variable (default \code{FALSE}).
#' @param verbose If \code{TRUE}, prints iteration logs.
#'
#' @return A \code{dimv_imputer} object containing fitted models, column means,
#' residual variances, convergence \code{diff_history}, and metadata.
#'
#' @examples
#' set.seed(1)
#' X <- data.frame(a = rnorm(50), b = rnorm(50), c = rnorm(50))
#' for (j in seq_along(X)) X[sample(1:nrow(X), 5), j] <- NA
#' imp <- dimv_train(X, lambda = 0.5)
#' print(imp)
#'
#' @export
#' @importFrom MASS ginv
#' @importFrom stats cor
dimv_train <- function(X, lambda = 1.0, maxit = 50, tol = 1e-4,
                       feature_select = FALSE, adaptive = FALSE, verbose = FALSE) {
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
  
  ridge_fit <- function(Xmat, y, lambda) {
    Z <- cbind(1, Xmat)
    D <- diag(c(0, rep(1, ncol(Xmat))))
    A <- crossprod(Z) + lambda * D
    b <- crossprod(Z, y)
    coef <- tryCatch(solve(A, b), error = function(e) MASS::ginv(A) %*% b)
    list(intercept = coef[1], beta = coef[-1])
  }
  
  prev <- X_imp
  models <- vector("list", length = p); names(models) <- coln
  diff_history <- numeric(0)
  it <- 0
  
  for (it in seq_len(maxit)) {
    for (j in seq_len(p)) {
      yvec <- prev[, j]
      
      # Feature selection (indices)
      if (feature_select) {
        if (exists("select_features_adaptive", mode = "function")) {
          sel <- try(select_features_adaptive(prev, target_var = j,
                                              min_features = 2, max_features = min(5, p - 1),
                                              verbose = FALSE), silent = TRUE)
          if (!inherits(sel, "try-error")) {
            if (is.character(sel$selected_features)) {
              feat_idx <- match(sel$selected_features, colnames(prev))
            } else {
              feat_idx <- as.integer(sel$selected_features)
            }
            feat_idx <- feat_idx[!is.na(feat_idx)]
          } else {
            feat_idx <- feature_select_fn(prev, j, max_features = min(5, p - 1))
          }
        } else {
          feat_idx <- feature_select_fn(prev, j, max_features = min(5, p - 1))
        }
        feat_idx <- intersect(feat_idx, setdiff(seq_len(p), j))
        if (length(feat_idx) == 0) feat_idx <- setdiff(seq_len(p), j)
        Xo <- prev[, feat_idx, drop = FALSE]
      } else {
        feat_idx <- setdiff(seq_len(p), j)
        Xo <- prev[, feat_idx, drop = FALSE]
      }
      
      lambda_j <- if (adaptive) adaptive_lambda_fn(Xo, base_lambda = lambda) else lambda
      
      fit <- ridge_fit(Xo, yvec, lambda_j)
      models[[j]] <- list(fit = fit, features = as.integer(feat_idx), lambda = lambda_j)
      
      na_idx <- which(is.na(as.matrix(X)[, j]))
      if (length(na_idx)) {
        preds <- fit$intercept + as.numeric(prev[na_idx, feat_idx, drop = FALSE] %*% fit$beta)
        prev[na_idx, j] <- preds
      }
    }
    
    if (any(nas)) {
      diff <- sqrt(mean((prev[nas] - X_imp[nas])^2))
      if (is.na(diff) || is.nan(diff)) diff <- 0
      diff_history <- c(diff_history, diff)
      if (verbose) cat("iter", it, "diff =", diff, "\n")
      X_imp <- prev
      if (diff < tol) break
    } else {
      if (verbose) cat("iter", it, ": no missing; stopping.\n")
      break
    }
  }
  
  resid_var <- numeric(p)
  for (j in seq_len(p)) {
    info <- models[[j]]
    fit <- info$fit
    Xo <- X_imp[, info$features, drop = FALSE]
    yvec <- X_imp[, j]
    predj <- fit$intercept + Xo %*% fit$beta
    resid_var[j] <- mean((yvec - predj)^2)
  }
  
  structure(list(models = models,
                 col_means = col_means,
                 resid_var = resid_var,
                 lambda = lambda,
                 adaptive = adaptive,
                 feature_select = feature_select,
                 iters = it,
                 diff_history = diff_history,
                 colnames = coln,
                 version = "0.1.2"),
            class = "dimv_imputer")
}

#' Impute New Data with a Fitted DIMV Imputer
#'
#' Applies a trained DIMV imputer to new (or the same) data. Iteratively imputes
#' missing values until convergence, using previously learned conditional models.
#'
#' @param imputer A \code{dimv_imputer} object returned by \code{dimv_train()}.
#' @param X_new A numeric data frame or matrix with missing values.
#' @param maxit Maximum number of update iterations (default \code{50}).
#' @param tol Convergence threshold (default \code{1e-4}).
#' @param verbose If \code{TRUE}, prints iteration logs.
#'
#' @return A fully imputed data frame.
#'
#' @export
dimv_impute_new <- function(imputer, X_new, maxit = 50, tol = 1e-4, verbose = FALSE) {
  if (!inherits(imputer, "dimv_imputer")) stop("imputer must be result of dimv_train()")
  X_new <- as.data.frame(lapply(X_new, as.numeric))
  X_new <- X_new[, imputer$colnames, drop = FALSE]
  
  X_imp <- as.matrix(X_new)
  nas <- is.na(X_imp)
  for (j in seq_along(imputer$col_means)) X_imp[is.na(X_imp[, j]), j] <- imputer$col_means[j]
  prev <- X_imp
  
  for (it in seq_len(maxit)) {
    for (j in seq_along(imputer$models)) {
      info <- imputer$models[[j]]
      fit <- info$fit
      feat_idx <- info$features
      Xo <- prev[, feat_idx, drop = FALSE]
      
      na_idx <- which(nas[, j])
      if (length(na_idx)) {
        preds <- fit$intercept + as.numeric(prev[na_idx, feat_idx, drop = FALSE] %*% fit$beta)
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
#' @param imputer A \code{dimv_imputer} model fitted by \code{dimv_train()}.
#' @param X_new Data containing missing values.
#' @param m Number of multiply-imputed datasets to generate (default \code{5}).
#' @param seed Optional random seed for reproducibility.
#'
#' @return A list of \code{m} fully-imputed data frames.
#'
#' @export
#' @importFrom stats rnorm
dimv_impute_multiple <- function(imputer, X_new, m = 5, seed = NULL) {
  if (!inherits(imputer, "dimv_imputer")) stop("imputer must be result of dimv_train()")
  if (!is.null(seed)) set.seed(seed)
  
  X_new <- as.data.frame(lapply(X_new, as.numeric))
  X_new <- X_new[, imputer$colnames, drop = FALSE]
  
  base_imp <- dimv_impute_new(imputer, X_new)
  
  res_list <- vector("list", m)
  for (k in seq_len(m)) {
    Xk <- as.matrix(base_imp)
    nas <- is.na(as.matrix(X_new))
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

#' Basic diagnostics for a trained DIMV imputer
#'
#' Provides a small set of diagnostics summarizing convergence, variance,
#' and a naive RMSE on previously-missing entries relative to column means.
#' Intended as a quick health check; not a replacement for downstream validation.
#'
#' @param imputer A \code{dimv_imputer} object.
#' @param X_original The original data (with \code{NA}).
#' @param X_imputed The imputed data returned by \code{dimv_impute_new()}.
#'
#' @return A list with \code{iterations}, \code{adaptive}, \code{feature_select},
#' \code{residual_variances}, \code{diff_history}, and \code{missing_rmse}.
#'
#' @export
dimv_diagnostics <- function(imputer, X_original, X_imputed) {
  nas <- is.na(as.matrix(X_original))
  rmse <- sqrt(mean((as.matrix(X_imputed)[nas] - imputer$col_means[col(X_imputed)[nas]])^2))
  list(
    iterations = imputer$iters,
    adaptive = imputer$adaptive,
    feature_select = imputer$feature_select,
    residual_variances = imputer$resid_var,
    diff_history = imputer$diff_history,
    missing_rmse = rmse
  )
}

#' @exportS3Method print dimv_imputer
print.dimv_imputer <- function(x, ...) {
  cat("R-native DIMV imputer (ridge-conditional)\n",
      "lambda =", x$lambda,
      "| adaptive =", x$adaptive,
      "| feature_select =", x$feature_select,
      "| iterations =", x$iters,
      "| version =", x$version, "\n")
}
