# ===========================================
# File: R/mean_helpers.R
# ===========================================

#' Fit Mean Imputer
#' 
#' Calculates column means for simple mean imputation
#' 
#' @param X_train Training data frame with missing values
#' @return List with column means
#' @keywords internal
fit_mean_imputer <- function(X_train) {
  list(means = sapply(X_train, function(x) mean(x, na.rm = TRUE)))
}

#' Apply Mean Imputer
#' 
#' Imputes missing values using pre-calculated means
#' 
#' @param imputer Result from fit_mean_imputer()
#' @param X Data frame to impute
#' @return Imputed data frame
#' @keywords internal
apply_mean_imputer <- function(imputer, X) {
  as.data.frame(
    mapply(
      function(col, m) { 
        col[is.na(col)] <- m
        col 
      }, 
      X, 
      imputer$means, 
      SIMPLIFY = FALSE
    )
  )
}