# ===========================================
# File: R/shap_parallel.R
# Robust Parallel SHAP Computation (Phase-1)
# ===========================================

#' Compute SHAP Values in Parallel
#'
#' Efficient parallelized SHAP value computation using fastshap + future.
#' Auto-fallback to sequential mode if parallel execution fails.
#'
#' @param model Trained model object (e.g., xgboost model)
#' @param X_train Training data frame used as background
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
  
  # Dependency validation
  for (pkg in c("xgboost", "fastshap", "future", "future.apply")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for compute_shap_parallel()", pkg))
    }
  }
  
  # Flexible prediction wrapper (fastshap version-agnostic)
  pred_wrapper <- function(..., newdata = NULL, object = NULL) {
    data <- if (!is.null(newdata)) newdata else object
    if (is.null(data)) stop("pred_wrapper: no data provided")
    as.numeric(predict(model, xgboost::xgb.DMatrix(as.matrix(data))))
  }
  
  # Split workload
  n <- nrow(X_test)
  splits <- split(seq_len(n), cut(seq_len(n), breaks = n_workers, labels = FALSE))
  
  # Plan management
  old_plan <- future::plan()                  # store current plan
  on.exit(future::plan(old_plan), add = TRUE) # restore after execution
  
  # Try parallel, fallback to sequential if needed
  try_parallel <- TRUE
  result <- NULL
  
  if (try_parallel) {
    try({
      future::plan(future::multisession, workers = n_workers)
      
      result <- future.apply::future_lapply(
        splits,
        function(idx) {
          requireNamespace("xgboost", quietly = TRUE)
          requireNamespace("fastshap", quietly = TRUE)
          
          args <- list(
            X = as.data.frame(X_train),
            pred_wrapper = pred_wrapper,
            nsim = nsim,
            newdata = as.data.frame(X_test[idx, , drop = FALSE]),
            adjust = TRUE
          )
          
          tryCatch(
            do.call(fastshap::explain, c(list(object = model), args)),
            error = function(e) do.call(fastshap::explain, args)
          )
        },
        future.seed = TRUE
      )
    }, silent = TRUE)
  }
  
  # If parallel failed => sequential fallback
  if (is.null(result)) {
    message("[SHAP] Parallel execution failed => running sequentially.")
    future::plan(future::sequential)
    
    result <- lapply(
      splits,
      function(idx) {
        args <- list(
          X = as.data.frame(X_train),
          pred_wrapper = pred_wrapper,
          nsim = nsim,
          newdata = as.data.frame(X_test[idx, , drop = FALSE]),
          adjust = TRUE
        )
        tryCatch(
          do.call(fastshap::explain, c(list(object = model), args)),
          error = function(e) do.call(fastshap::explain, args)
        )
      }
    )
  }
  
  do.call(rbind, result)
}
