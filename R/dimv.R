# ===========================================
# File: R/dimv.R
# ===========================================
# R-native DIMV (ridge-conditional) with multiple-imputation helper

dimv_train <- function(X, lambda = 1.0, maxit = 50, tol = 1e-4, verbose = FALSE) {
  X <- as.data.frame(lapply(X, as.numeric))
  n <- nrow(X); p <- ncol(X)
  coln <- colnames(X)
  if (is.null(coln)) coln <- paste0("V", seq_len(p))
  colnames(X) <- coln
  
  # initial imputation
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
  it <- 0
  
  for (it in seq_len(maxit)) {
    for (j in seq_len(p)) {
      yvec <- prev[, j]
      Xo <- prev[, -j, drop = FALSE]
      fit <- ridge_fit(Xo, yvec, lambda)
      models[[j]] <- fit
      na_idx <- which(is.na(as.matrix(X)[, j]))
      if (length(na_idx)) {
        preds <- fit$intercept + as.numeric(Xo[na_idx, , drop = FALSE] %*% fit$beta)
        prev[na_idx, j] <- preds
      }
    }
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
  
  resid_var <- numeric(p)
  for (j in seq_len(p)) {
    yvec <- X_imp[, j]
    Xo <- X_imp[, -j, drop = FALSE]
    predj <- models[[j]]$intercept + Xo %*% models[[j]]$beta
    resid_var[j] <- mean((yvec - predj)^2)
  }
  
  structure(list(models = models,
                 col_means = col_means,
                 resid_var = resid_var,
                 lambda = lambda,
                 iters = it,
                 colnames = coln),
            class = "dimv_imputer")
}

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

# -----------------------------
# Multiple-imputation helper
# -----------------------------
# Generate m imputed datasets using DIMV by adding conditional noise from resid_var.
dimv_impute_multiple <- function(imputer, X_new, m = 5, seed = NULL) {
  if (!inherits(imputer, "dimv_imputer")) stop("imputer must be result of dimv_train()")
  if (!is.null(seed)) set.seed(seed)
  X_new <- as.data.frame(lapply(X_new, as.numeric))
  n <- nrow(X_new); p <- ncol(X_new)
  coln <- imputer$colnames
  X_new <- X_new[, coln, drop = FALSE]
  base_imp <- dimv_impute_new(imputer, X_new)  # deterministic conditional expectation
  res_list <- vector("list", m)
  for (k in seq_len(m)) {
    Xk <- as.matrix(base_imp)
    # add Gaussian noise only to entries that were originally NA in X_new
    nas <- is.na(as.matrix(X_new))
    for (j in seq_len(p)) {
      if (any(nas[, j])) {
        sigma <- sqrt(max(0, imputer$resid_var[j]))
        Xk[nas[, j], j] <- Xk[nas[, j], j] + rnorm(sum(nas[, j]), mean = 0, sd = sigma)
      }
    }
    res_list[[k]] <- as.data.frame(Xk)
  }
  res_list
}

#' @export
print.dimv_imputer <- function(x, ...) {
  cat("R-native DIMV imputer (ridge-conditional)\n",
      "lambda =", x$lambda, "| iterations =", x$iters, "\n")
}
