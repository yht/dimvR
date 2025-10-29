# ===========================================
# File: R/shap_parallel.R
# ===========================================

compute_shap_parallel <- function(model, X_train, X_test, nsim = 1000, n_workers = 4) {
  # --- pastikan dependency tersedia ---
  requireNamespace("xgboost", quietly = TRUE)
  requireNamespace("fastshap", quietly = TRUE)
  requireNamespace("future", quietly = TRUE)
  requireNamespace("future.apply", quietly = TRUE)
  
  # --- pembungkus prediksi fleksibel ---
  # kompatibel untuk fastshap versi lama (object,newdata) dan baru (newdata)
  pred_wrapper <- function(..., newdata = NULL, object = NULL) {
    if (!is.null(newdata)) {
      return(as.numeric(predict(model, xgboost::xgb.DMatrix(as.matrix(newdata)))))
    } else if (!is.null(object)) {
      return(as.numeric(predict(model, xgboost::xgb.DMatrix(as.matrix(object)))))
    } else {
      stop("pred_wrapper: tidak ada data diberikan")
    }
  }
  
  # --- parallel setup ---
  n <- nrow(X_test)
  splits <- split(seq_len(n), cut(seq_len(n), breaks = n_workers, labels = FALSE))
  
  future::plan(future::multisession, workers = n_workers)
  
  parts <- future.apply::future_lapply(
    splits,
    function(idx) {
      # pastikan dependency juga aktif di tiap worker
      requireNamespace("xgboost", quietly = TRUE)
      requireNamespace("fastshap", quietly = TRUE)
      
      args <- list(
        X = as.data.frame(X_train),
        pred_wrapper = pred_wrapper,
        nsim = nsim,
        newdata = as.data.frame(X_test[idx, , drop = FALSE]),
        adjust = TRUE
      )
      
      # fallback kompatibilitas fastshap lama
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
