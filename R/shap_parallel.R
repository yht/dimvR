# ===========================================
# File: R/shap_parallel.R
# ===========================================

#' Compute SHAP Values in Parallel
#' 
#' Calculates SHAP (SHapley Additive exPlanations) values for test data
#' using parallel processing for efficiency.
#' 
#' @param model Trained model object (e.g., xgboost model)
#' @param X_train Training data frame used as background data
#' @param X_test Test data frame to explain
#' @param nsim Number of Monte Carlo simulations for SHAP
#' @param n_workers Number of parallel workers
#' @return Matrix of SHAP values
#' @keywords internal
#' @importFrom xgboost xgb.DMatrix
#' @importFrom fastshap explain
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom stats predict
compute_shap_parallel <- function(model, X_train, X_test, nsim = 1000, n_workers = 4) {
  # Check dependencies
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required for compute_shap_parallel()")
  }
  if (!requireNamespace("fastshap", quietly = TRUE)) {
    stop("Package 'fastshap' is required for compute_shap_parallel()")
  }
  if (!requireNamespace("future", quietly = TRUE)) {
    stop("Package 'future' is required for compute_shap_parallel()")
  }
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Package 'future.apply' is required for compute_shap_parallel()")
  }
  
  # Flexible prediction wrapper
  # Compatible with both old (object, newdata) and new (newdata) fastshap versions
  pred_wrapper <- function(..., newdata = NULL, object = NULL) {
    if (!is.null(newdata)) {
      return(as.numeric(predict(model, xgboost::xgb.DMatrix(as.matrix(newdata)))))
    } else if (!is.null(object)) {
      return(as.numeric(predict(model, xgboost::xgb.DMatrix(as.matrix(object)))))
    } else {
      stop("pred_wrapper: no data provided")
    }
  }
  
  # Parallel setup
  n <- nrow(X_test)
  splits <- split(seq_len(n), cut(seq_len(n), breaks = n_workers, labels = FALSE))
  
  future::plan(future::multisession, workers = n_workers)
  
  parts <- future.apply::future_lapply(
    splits,
    function(idx) {
      # Ensure dependencies are loaded in each worker
      requireNamespace("xgboost", quietly = TRUE)
      requireNamespace("fastshap", quietly = TRUE)
      
      args <- list(
        X = as.data.frame(X_train),
        pred_wrapper = pred_wrapper,
        nsim = nsim,
        newdata = as.data.frame(X_test[idx, , drop = FALSE]),
        adjust = TRUE
      )
      
      # Fallback for old fastshap versions
      tryCatch({
        do.call(fastshap::explain, c(list(object = model), args))
      }, error = function(e) {
        do.call(fastshap::explain, args)
      })
    },
    future.seed = TRUE
  )
  
  future::plan(future::sequential)
  do.call(rbind, parts)
}