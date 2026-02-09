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
#' @param seed Optional RNG seed for reproducible SHAP estimates
#' @return Matrix of SHAP values
#' @keywords internal
#' @importFrom xgboost xgb.DMatrix
#' @importFrom fastshap explain
#' @importFrom future plan multisession sequential
#' @importFrom future.apply future_lapply
#' @importFrom stats predict
compute_shap_parallel <- function(model, X_train, X_test, nsim = 1000, n_workers = 4, seed = NULL) {
  
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
    mat <- as.matrix(data)
    # Newer xgboost interfaces prefer matrix input for class 'xgboost';
    # keep DMatrix fallback for older versions/backends.
    pred <- tryCatch(
      stats::predict(model, mat),
      error = function(e) stats::predict(model, xgboost::xgb.DMatrix(mat))
    )
    as.numeric(pred)
  }
  
  # Split workload (ensure consistent chunking across sequential/parallel)
  n <- nrow(X_test)
  if (n <= 1) {
    splits <- list(seq_len(n))
  } else {
    split_count <- if (is.null(n_workers) || !is.numeric(n_workers) || n_workers <= 1) {
      2
    } else {
      n_workers
    }
    split_count <- min(split_count, n)
    splits <- split(seq_len(n), cut(seq_len(n), breaks = split_count, labels = FALSE))
  }
  
  # Plan management
  old_plan <- future::plan()                  # store current plan
  on.exit(future::plan(old_plan), add = TRUE) # restore after execution
  
  # Try parallel, fallback to sequential if needed
  try_parallel <- is.numeric(n_workers) && n_workers > 1
  result <- NULL
  
  seed_use <- if (is.null(seed)) 123 else seed
  
  if (try_parallel) {
    try({
      future::plan(future::multisession, workers = n_workers)
      
      split_seeds <- seed_use + seq_along(splits)
      result <- future.apply::future_lapply(
        seq_along(splits),
        function(i) {
          set.seed(split_seeds[i])
          idx <- splits[[i]]
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
        future.seed = seed_use
      )
    }, silent = TRUE)
  }
  
  # If parallel failed => sequential fallback
  if (is.null(result)) {
    message("[SHAP] Parallel execution failed => running sequentially.")
    future::plan(future::sequential)
    
    split_seeds <- seed_use + seq_along(splits)
    result <- future.apply::future_lapply(
      seq_along(splits),
      function(i) {
        set.seed(split_seeds[i])
        idx <- splits[[i]]
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
      future.seed = seed_use
    )
  }
  
  do.call(rbind, result)
}
