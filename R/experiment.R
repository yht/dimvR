# R/experiment.R
# run_full_pipeline with multiple-imputation pooling (Rubin)

# helper: compute within-imputation variance for MSE estimator
# mse = mean(e_i^2); Var(mse) ≈ Var(e_i^2)/n
within_var_mse <- function(y_true, y_hat) {
  errs2 <- (as.numeric(y_true) - as.numeric(y_hat))^2
  n <- length(errs2)
  if (n <= 1) return(0)
  var(errs2, na.rm = TRUE) / n
}

# within-var for SHAP MSE (treat all entries equally)
within_var_shap_mse <- function(shap_mat, ref_shap) {
  # flatten elementwise squared errors
  se <- as.numeric((as.matrix(shap_mat) - as.matrix(ref_shap))^2)
  N <- length(se)
  if (N <= 1) return(0)
  var(se, na.rm = TRUE) / N
}

# Rubin pooling: Q_m vector, U_m vector
rubin_pool <- function(Qm, Um) {
  m <- length(Qm)
  Q_bar <- mean(Qm)
  B <- if (m > 1) var(Qm) else 0
  U_bar <- mean(Um)
  T_var <- U_bar + (1 + 1/m) * B
  se_total <- sqrt(T_var)
  list(Q_bar = Q_bar, se = se_total, between = B, within = U_bar, m = m)
}

# simulate missing helper (same as previous)
simulate_missing <- function(X, rate = 0.2, mechanism = "MCAR") {
  X <- as.data.frame(X)
  n <- nrow(X); p <- ncol(X)
  M <- matrix(FALSE, n, p)
  if (mechanism == "MCAR") {
    M <- matrix(runif(n * p) < rate, n, p)
  } else if (mechanism == "MAR") {
    prob <- plogis(scale(X[, 1], center = TRUE, scale = TRUE))
    for (j in seq_len(p)) M[, j] <- runif(n) < (rate * prob)
  } else if (mechanism == "MNAR") {
    for (j in seq_len(p)) {
      prob <- plogis(scale(X[, j], center = TRUE, scale = TRUE))
      M[, j] <- runif(n) < (rate * prob)
    }
  } else stop("Unknown mechanism: ", mechanism)
  X_missing <- X; X_missing[M] <- NA
  X_missing
}

run_full_pipeline <- function(X, y,
                              task = "regression",
                              missing_rates = c(0.2),
                              mechanisms = c("MCAR"),
                              imputers = c("dimv_R", "mean", "mice"),
                              nsim = 1000,
                              repeats = 1,
                              workers = 2,
                              m_imp = 5,      # number of imputations for pooling
                              seed = 12345) {
  
  requireNamespace("xgboost", quietly = TRUE)
  requireNamespace("future", quietly = TRUE)
  requireNamespace("fastshap", quietly = TRUE)
  requireNamespace("mice", quietly = TRUE)
  
  set.seed(seed)
  results_list <- list()
  run_id <- 1
  
  for (rate in missing_rates) {
    for (mech in mechanisms) {
      cat("=== rate:", rate, "mechanism:", mech, "===\n")
      X_miss_full <- simulate_missing(X, rate = rate, mechanism = mech)
      X_miss_full <- as.data.frame(lapply(X_miss_full, as.numeric))
      
      # split by rows (use original y without NA)
      tr_idx <- caret::createDataPartition(y, p = 0.8, list = FALSE)
      X_train_orig <- X[tr_idx, , drop = FALSE]
      X_test_orig  <- X[-tr_idx, , drop = FALSE]
      y_train_orig <- y[tr_idx]; y_test_orig <- y[-tr_idx]
      
      # reference (full-data) model and SHAP
      ref_mod <- xgboost::xgboost(data = as.matrix(X_train_orig), label = y_train_orig,
                                  nrounds = 200, objective = "reg:squarederror", verbose = 0)
      ref_shap <- compute_shap_parallel(ref_mod, X_train_orig, X_test_orig, nsim = nsim, n_workers = workers)
      
      # now create missing train/test from X_miss_full split
      X_train_miss <- X_miss_full[tr_idx, , drop = FALSE]
      X_test_miss  <- X_miss_full[-tr_idx, , drop = FALSE]
      
      for (imp in imputers) {
        cat("Processing imputer:", imp, "...\n")
        # storage for imputation-specific estimates
        Q_m_pred <- numeric(m_imp)
        U_m_pred <- numeric(m_imp)
        Q_m_shap <- numeric(m_imp)
        U_m_shap <- numeric(m_imp)
        
        if (imp == "mice") {
          # fit mice on train, get m completed trainings, impute test via model-based predict per completed dataset
          mice_obj <- mice::mice(X_train_miss, m = m_imp, maxit = 5, printFlag = FALSE)
          completed_trains <- lapply(1:m_imp, function(k) mice::complete(mice_obj, k))
          # create pooled test imputations: here we will produce m imputed tests via simple predictive models per completed train
          imp_tests <- vector("list", m_imp)
          for (k in seq_len(m_imp)) {
            trk <- completed_trains[[k]]
            # for test imputation we use fitted RF per-variable as in earlier code
            imp_tests[[k]] <- apply_mice_pool_imputer(list(miceobj = mice_obj, m = m_imp), X_train_miss, X_test_miss)
            # training data for model training
            Xtr_k <- trk
            # train model and eval per-imputation below
            mod_k <- xgboost::xgboost(data = as.matrix(Xtr_k), label = y_train_orig,
                                      nrounds = 200, objective = "reg:squarederror", verbose = 0)
            pred_k <- predict(mod_k, xgboost::xgb.DMatrix(as.matrix(imp_tests[[k]])))
            Q_m_pred[k] <- mean((pred_k - y_test_orig)^2)
            U_m_pred[k] <- within_var_mse(y_test_orig, pred_k)
            # SHAP
            shap_k <- compute_shap_parallel(mod_k, Xtr_k, imp_tests[[k]], nsim = nsim, n_workers = workers)
            Q_m_shap[k] <- mean((as.numeric(shap_k) - as.numeric(ref_shap))^2)
            U_m_shap[k] <- within_var_shap_mse(shap_k, ref_shap)
          }
        } else if (imp == "dimv_R") {
          # fit single dimv on train (conditional expectation) then generate m imputations by adding noise
          imputer <- dimv_train(X_train_miss)
          imp_tests <- dimv_impute_multiple(imputer, X_test_miss, m = m_imp, seed = seed)
          imp_trains <- dimv_impute_multiple(imputer, X_train_miss, m = m_imp, seed = seed + 1)
          for (k in seq_len(m_imp)) {
            Xtr_k <- imp_trains[[k]]
            Xte_k <- imp_tests[[k]]
            mod_k <- xgboost::xgboost(data = as.matrix(Xtr_k), label = y_train_orig,
                                      nrounds = 200, objective = "reg:squarederror", verbose = 0)
            pred_k <- predict(mod_k, xgboost::xgb.DMatrix(as.matrix(Xte_k)))
            Q_m_pred[k] <- mean((pred_k - y_test_orig)^2)
            U_m_pred[k] <- within_var_mse(y_test_orig, pred_k)
            shap_k <- compute_shap_parallel(mod_k, Xtr_k, Xte_k, nsim = nsim, n_workers = workers)
            Q_m_shap[k] <- mean((as.numeric(shap_k) - as.numeric(ref_shap))^2)
            U_m_shap[k] <- within_var_shap_mse(shap_k, ref_shap)
          }
        } else if (imp == "mean") {
          # single deterministic imputation (no pooling) -> treat as m=1
          imp_mean <- fit_mean_imputer(X_train_miss)
          Xtr_imp <- apply_mean_imputer(imp_mean, X_train_miss)
          Xte_imp <- apply_mean_imputer(imp_mean, X_test_miss)
          mod <- xgboost::xgboost(data = as.matrix(Xtr_imp), label = y_train_orig,
                                  nrounds = 200, objective = "reg:squarederror", verbose = 0)
          pred <- predict(mod, xgboost::xgb.DMatrix(as.matrix(Xte_imp)))
          Q_m_pred[1] <- mean((pred - y_test_orig)^2); U_m_pred[1] <- within_var_mse(y_test_orig, pred)
          shap_m <- compute_shap_parallel(mod, Xtr_imp, Xte_imp, nsim = nsim, n_workers = workers)
          Q_m_shap[1] <- mean((as.numeric(shap_m) - as.numeric(ref_shap))^2); U_m_shap[1] <- within_var_shap_mse(shap_m, ref_shap)
          if (m_imp > 1) {
            Q_m_pred <- rep(Q_m_pred[1], m_imp); U_m_pred <- rep(U_m_pred[1], m_imp)
            Q_m_shap <- rep(Q_m_shap[1], m_imp); U_m_shap <- rep(U_m_shap[1], m_imp)
          }
        } else {
          stop("Unknown imputer: ", imp)
        }
        
        # apply Rubin pooling
        pool_pred <- rubin_pool(Q_m_pred, U_m_pred)
        pool_shap <- rubin_pool(Q_m_shap, U_m_shap)
        
        results_list[[length(results_list) + 1]] <- data.frame(
          method = imp,
          mechanism = mech,
          rate = rate,
          mse_pred = pool_pred$Q_bar,
          mse_pred_se = pool_pred$se,
          mse_shap = pool_shap$Q_bar,
          mse_shap_se = pool_shap$se,
          m = pool_pred$m,
          stringsAsFactors = FALSE
        )
        
        run_id <- run_id + 1
      } # end imputer loop
    } # end mechanism loop
  } # end rate loop
  
  out <- do.call(rbind, results_list)
  rownames(out) <- NULL
  out
}
