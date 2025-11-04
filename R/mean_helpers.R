# ===========================================
# File: R/mean_helpers.R
# Mean Imputation (Baseline Comparator)
# Fully deterministic, zero-dependency
# ===========================================

#' Fit Mean Imputer (Deterministic Baseline)
#'
#' Computes per-column means on TRAINING data only,
#' used as a baseline method for benchmarking against DIMV/MICE.
#'
#' @param X_train Training data frame (may contain NA).
#' @return List containing column means and column names.
#' @keywords internal
fit_mean_imputer <- function(X_train) {
  if (!is.data.frame(X_train)) stop("X_train must be a data.frame")
  
  # Convert to numeric where possible (warn if not)
  Xn <- suppressWarnings(lapply(X_train, as.numeric))
  non_num <- which(vapply(Xn, function(x) all(is.na(x)), logical(1)))
  if (length(non_num) > 0) {
    warning("Non-numeric columns detected in mean imputation: ",
            paste(names(X_train)[non_num], collapse = ", "),
            ". Returning NA means for those columns.")
  }
  
  list(
    means = vapply(Xn, function(x) mean(x, na.rm = TRUE), numeric(1)),
    cols  = colnames(X_train)
  )
}

#' Apply Mean Imputer on New Data
#'
#' Imputes missing values using TRAIN-derived column means only.
#'
#' @param imputer Object returned by `fit_mean_imputer()`.
#' @param X Data frame to impute.
#' @return Imputed data frame (same column order as input).
#' @keywords internal
apply_mean_imputer <- function(imputer, X) {
  if (!is.list(imputer) || is.null(imputer$means)) {
    stop("Invalid imputer: must come from fit_mean_imputer()")
  }
  
  X_imp <- as.data.frame(X)
  
  # Align columns
  common <- intersect(names(imputer$means), names(X_imp))
  for (col in common) {
    nas <- is.na(X_imp[[col]])
    if (any(nas)) X_imp[[col]][nas] <- imputer$means[[col]]
  }
  
  X_imp
}
