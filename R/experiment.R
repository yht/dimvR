# ===========================================
# File: R/experiment.R
# run_full_pipeline with multiple-imputation pooling (Rubin)
# Requires: mice_backend_train(), mice_backend_apply() (in R/mice_backend.R)
# ===========================================

# ---------- Internal utilities (variance, pooling, miss mech) ----------

#' Compute Within-Imputation Variance for MSE of predictions
#' @keywords internal
within_var_mse <- function(y_true, y_hat) {
  e2 <- (as.numeric(y_true) - as.numeric(y_hat))^2
  n <- length(e2)
  if (n <= 1) return(0)
  stats::var(e2, na.rm = TRUE) / n
}

#' Within-variance for SHAP MSE (elementwise over matrix)
#' @keywords internal
within_var_shap_mse <- function(shap_mat, ref_shap) {
  A <- as.matrix(shap_mat)
  B <- as.matrix(ref_shap)
  if (!all(dim(A) == dim(B))) stop("within_var_shap_mse: SHAP matrices must have same shape.")
  se <- as.numeric((A - B)^2)
  N <- length(se)
  if (N <= 1) return(0)
  stats::var(se, na.rm = TRUE) / N
}

#' Rubin's Pooling Rules
#' @keywords internal
rubin_pool <- function(Qm, Um) {
  m <- length(Qm)
  Q_bar <- mean(Qm)
  B <- if (m > 1) stats::var(Qm) else 0
  U_bar <- mean(Um)
  T_var <- U_bar + (1 + 1/m) * B
  list(Q_bar = Q_bar,
       se = sqrt(T_var),
       between = B,
       within = U_bar,
       m = m)
}

#' Simulate Missing Data (MCAR/MAR/MNAR) safely
#' @keywords internal
simulate_missing <- function(X, rate = 0.2, mechanism = "MCAR") {
  X <- as.data.frame(lapply(X, as.numeric))
  n <- nrow(X); p <- ncol(X)
  M <- matrix(FALSE, n, p)
  clamp <- function(z) pmin(pmax(z, 0), 1)
  if (mechanism == "MCAR") {
    M <- matrix(stats::runif(n * p) < rate, n, p)
  } else if (mechanism == "MAR") {
    aux <- X[[1]]
    if (all(is.na(aux))) aux <- rowMeans(X, na.rm = TRUE)
    aux[is.na(aux)] <- mean(aux, na.rm = TRUE)
    prob_row <- clamp(stats::plogis(scale(aux)) * rate)
    for (j in seq_len(p)) M[, j] <- stats::runif(n) < prob_row
  } else if (mechanism == "MNAR") {
    for (j in seq_len(p)) {
      v <- X[[j]]
      v_tmp <- v
      v_tmp[is.na(v_tmp)] <- mean(v_tmp, na.rm = TRUE)
      pj <- clamp(stats::plogis(scale(v_tmp)) * rate)
      M[, j] <- stats::runif(n) < pj
    }
  } else {
    stop("Unknown mechanism: ", mechanism)
  }
  X_miss <- X
  X_miss[M] <- NA
  X_miss
}

# ---------- Package presence checks ----------

.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
}

# ---------- Main pipeline ----------

#' Run Full DIMV Pipeline with Multiple Imputation (Rubin pooling)
#'
#' Executes an end-to-end experiment:
#' 1) Induces missingness (MCAR/MAR/MNAR) at given rates
#' 2) Trains reference model on complete TRAIN and computes reference SHAP on TEST
#' 3) Imputes TRAIN/TEST via selected imputers (train-only where applicable)
#' 4) Trains downstream model per imputation, evaluates MSE & SHAP MSE
#' 5) Pools estimates with Rubin's rules
#'
#' @param X Data frame/matrix of predictors (numeric).
#' @param y Numeric response.
#' @param task "regression" (placeholder for future extensions).
#' @param missing_rates Vector of missing rates (e.g., c(0.2, 0.4)).
#' @param mechanisms Vector of mechanisms ("MCAR","MAR","MNAR").
#' @param imputers Vector of imputer keys: c("dimv_R","mean","mice").
#' @param nsim Monte Carlo simulations for SHAP (fastshap).
#' @param repeats Repeats per configuration (currently unused).
#' @param workers Number of parallel workers for SHAP.
#' @param m_imp Number of imputations per method (Rubin pooling).
#' @param seed RNG seed.
#'
#' @return Data frame with pooled metrics per (method, mechanism, rate).
#' @export
#' @importFrom stats plogis var runif predict
#' @importFrom xgboost xgboost xgb.DMatrix
run_full_pipeline <- function(X, y,
                              task = "regression",
                              missing_rates = c(0.2),
                              mechanisms = c("MCAR"),
                              imputers = c("dimv_R", "mean", "mice"),
                              nsim = 1000,
                              repeats = 1,
                              workers = 2,
                              m_imp = 5,
                              seed = 12345) {
  
  .require_pkg("xgboost")
  .require_pkg("fastshap")
  # MICE only required if requested
  if ("mice" %in% imputers) .require_pkg("mice")
  
  set.seed(seed)
  X <- as.data.frame(lapply(X, as.numeric))
  y <- as.numeric(y)
  
  results <- list()
  
  # ---------- Data split (base R) ----------
  n <- nrow(X)
  idx <- sample.int(n, size = floor(0.8 * n))
  X_train_orig <- X[idx, , drop = FALSE]
  X_test_orig  <- X[-idx, , drop = FALSE]
  y_train_orig <- y[idx]
  y_test_orig  <- y[-idx]
  
  # ---------- Reference model & SHAP on fully observed data ----------
  ref_mod <- xgboost::xgboost(
    data = as.matrix(X_train_orig),
    label = y_train_orig,
    nrounds = 200,
    objective = "reg:squarederror",
    verbose = 0
  )
  ref_shap <- compute_shap_parallel(ref_mod, X_train_orig, X_test_orig,
                                    nsim = nsim, n_workers = workers)
  
  for (rate in missing_rates) {
    for (mech in mechanisms) {
      
      cat("=== rate:", rate, "| mechanism:", mech, "===\n")
      X_miss_full <- simulate_missing(X, rate = rate, mechanism = mech)
      
      X_train_miss <- X_miss_full[idx, , drop = FALSE]
      X_test_miss  <- X_miss_full[-idx, , drop = FALSE]
      
      for (imp in imputers) {
        cat(" => Imputer:", imp, "\n")
        
        Q_m_pred <- numeric(m_imp)  # pooled estimands (per draw)
        U_m_pred <- numeric(m_imp)  # within-imputation variances
        Q_m_shap <- numeric(m_imp)
        U_m_shap <- numeric(m_imp)
        
        if (imp == "dimv_R") {
          # Train DIMV on TRAIN only; sample m imputations for TRAIN & TEST
          imp_obj    <- dimv_train(X_train_miss)
          imp_trains <- dimv_impute_multiple(imp_obj, X_train_miss, m = m_imp, seed = seed + 1)
          imp_tests  <- dimv_impute_multiple(imp_obj, X_test_miss,  m = m_imp, seed = seed + 2)
          
          for (k in seq_len(m_imp)) {
            Xtr_k <- imp_trains[[k]]
            Xte_k <- imp_tests[[k]]
            
            mod_k <- xgboost::xgboost(
              data = as.matrix(Xtr_k),
              label = y_train_orig,
              nrounds = 200,
              objective = "reg:squarederror",
              verbose = 0
            )
            pred_k <- stats::predict(mod_k, xgboost::xgb.DMatrix(as.matrix(Xte_k)))
            Q_m_pred[k] <- mean((pred_k - y_test_orig)^2)
            U_m_pred[k] <- within_var_mse(y_test_orig, pred_k)
            
            shap_k <- compute_shap_parallel(mod_k, Xtr_k, Xte_k, nsim = nsim, n_workers = workers)
            Q_m_shap[k] <- mean((as.matrix(shap_k) - as.matrix(ref_shap))^2)
            U_m_shap[k] <- within_var_shap_mse(shap_k, ref_shap)
          }
          
        } else if (imp == "mean") {
          # Train mean-imputer on TRAIN only, apply to TRAIN/TEST
          imp_mean <- fit_mean_imputer(X_train_miss)
          Xtr_imp  <- apply_mean_imputer(imp_mean, X_train_miss)
          Xte_imp  <- apply_mean_imputer(imp_mean, X_test_miss)
          
          mod <- xgboost::xgboost(
            data = as.matrix(Xtr_imp),
            label = y_train_orig,
            nrounds = 200,
            objective = "reg:squarederror",
            verbose = 0
          )
          pred <- stats::predict(mod, xgboost::xgb.DMatrix(as.matrix(Xte_imp)))
          mse_pred <- mean((pred - y_test_orig)^2)
          Q_m_pred[] <- mse_pred
          U_m_pred[] <- within_var_mse(y_test_orig, pred)
          
          shap_m <- compute_shap_parallel(mod, Xtr_imp, Xte_imp, nsim = nsim, n_workers = workers)
          mse_shap <- mean((as.matrix(shap_m) - as.matrix(ref_shap))^2)
          Q_m_shap[] <- mse_shap
          U_m_shap[] <- within_var_shap_mse(shap_m, ref_shap)
          
        } else if (imp == "mice") {
          # ---- Backend abstraction (train-only semantics) ----
          # mice_backend_train(): fit chained equations on TRAIN only
          # mice_backend_apply():  apply learned structure to new TEST data
          if (!exists("mice_backend_train", mode = "function") ||
              !exists("mice_backend_apply", mode = "function")) {
            stop("mice backend not found. Please add R/mice_backend.R with mice_backend_train/apply.")
          }
          
          mice_fit   <- mice_backend_train(X_train_miss, m = m_imp, maxit = 5, seed = seed + 10)
          imp_trains <- mice_fit$imputations
          imp_tests  <- lapply(seq_len(m_imp), function(k) {
            mice_backend_apply(mice_fit, X_test_miss, k = k)
          })
          
          for (k in seq_len(m_imp)) {
            Xtr_k <- imp_trains[[k]]
            Xte_k <- imp_tests[[k]]
            
            mod_k <- xgboost::xgboost(
              data = as.matrix(Xtr_k),
              label = y_train_orig,
              nrounds = 200,
              objective = "reg:squarederror",
              verbose = 0
            )
            pred_k <- stats::predict(mod_k, xgboost::xgb.DMatrix(as.matrix(Xte_k)))
            Q_m_pred[k] <- mean((pred_k - y_test_orig)^2)
            U_m_pred[k] <- within_var_mse(y_test_orig, pred_k)
            
            shap_k <- compute_shap_parallel(mod_k, Xtr_k, Xte_k, nsim = nsim, n_workers = workers)
            Q_m_shap[k] <- mean((as.matrix(shap_k) - as.matrix(ref_shap))^2)
            U_m_shap[k] <- within_var_shap_mse(shap_k, ref_shap)
          }
          
        } else {
          stop("Unknown imputer: ", imp)
        }
        
        # Rubin pooling
        pool_pred <- rubin_pool(Q_m_pred, U_m_pred)
        pool_shap <- rubin_pool(Q_m_shap, U_m_shap)
        
        results[[length(results) + 1]] <- data.frame(
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
      }
    }
  }
  
  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}
