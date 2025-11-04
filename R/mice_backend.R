# ===========================================
# File: R/mice_backend.R
# MICE Backend Abstraction (M2 Architecture)
#
# Provides TRAIN-only MICE fitting and APPLY-to-TEST functionality
# without leaking test-set information.
#
# This is an EXPERIMENTAL backend.
# Export only when API stabilized.
# ===========================================

#' Train MICE Backend (TRAIN-only)
#'
#' Fits MICE on TRAIN data only and returns `m` completed TRAIN datasets
#' plus the fitted MICE object, to allow application on TEST data later.
#'
#' @param X_train Training data.frame with missing values
#' @param m Number of imputations
#' @param maxit Number of MICE iterations
#' @param seed RNG seed
#' @param ... Additional args passed to `mice::mice()`
#'
#' @return A list with:
#'   - imputations: list of length `m` containing completed TRAIN datasets
#'   - model: fitted mice object
#'   - m: number of imputations
#'
#' @keywords internal
mice_backend_train <- function(X_train, m = 5, maxit = 5, seed = NULL, ...) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required for mice_backend_train()")
  }
  if (!is.null(seed)) set.seed(seed)
  
  X_train <- as.data.frame(X_train)
  
  mice_fit <- mice::mice(
    data   = X_train,
    m      = m,
    maxit  = maxit,
    printFlag = FALSE,
    ...
  )
  
  imputations <- lapply(seq_len(m), function(k) mice::complete(mice_fit, k))
  names(imputations) <- paste0("imp", seq_len(m))
  
  list(
    imputations = imputations,
    model = mice_fit,
    m = m
  )
}


#' Apply MICE Backend to TEST Data
#'
#' Applies a TRAIN-fitted MICE model to TEST data in a compatible way.
#' TEST data is imputed by training a model for each variable using the
#' corresponding completed TRAIN dataset. This keeps TRAIN-only semantics.
#'
#' @param mice_fit_list List returned by mice_backend_train()
#' @param X_test Test data.frame with missing values
#' @param k Which imputation index to apply (1..m)
#'
#' @return Imputed TEST data.frame
#' @keywords internal
mice_backend_apply <- function(mice_fit_list, X_test, k = 1) {
  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required for mice_backend_apply()")
  }
  
  X_test_imp <- as.data.frame(X_test)
  train_k <- mice::complete(mice_fit_list$model, k)
  
  for (col in names(X_test_imp)) {
    miss_idx <- which(is.na(X_test_imp[[col]]))
    if (!length(miss_idx)) next
    
    preds <- setdiff(names(train_k), col)
    train_df <- cbind(y = train_k[[col]], train_k[, preds, drop = FALSE])
    train_df <- train_df[stats::complete.cases(train_df), , drop = FALSE]
    
    if (nrow(train_df) < 10) {
      mu <- mean(train_df$y, na.rm = TRUE)
      X_test_imp[miss_idx, col] <- mu
    } else {
      # simple linear model for speed & CRAN safety
      fit <- try(stats::lm(y ~ ., data = train_df), silent = TRUE)
      if (inherits(fit, "try-error")) {
        mu <- mean(train_df$y, na.rm = TRUE)
        X_test_imp[miss_idx, col] <- mu
      } else {
        preds_test <- X_test_imp[miss_idx, preds, drop = FALSE]
        preds_test[is.na(preds_test)] <- 0
        X_test_imp[miss_idx, col] <- stats::predict(fit, preds_test)
      }
    }
  }
  
  X_test_imp
}
