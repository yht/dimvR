# ===========================================
# File: R/mice_helpers.R
# ===========================================
# Apply fitted mice object to new data
# We extract the predictive mean matching models from the fitted mice object
# and use them to fill missing values in newdata.

apply_mice_pool_imputer <- function(miceobj_list, X_train, X_test) {
  if (!requireNamespace("mice", quietly = TRUE))
    stop("Package 'mice' must be installed to use apply_mice_pool_imputer().")
  
  miceobj <- miceobj_list$miceobj
  m <- miceobj_list$m
  
  X_test_imp <- as.data.frame(X_test)
  X_train_comp <- mice::complete(miceobj, 1)
  
  # For each variable with NA in test, train a simple regression on train to predict it
  for (j in names(X_test_imp)) {
    nas <- is.na(X_test_imp[[j]])
    if (any(nas)) {
      # choose predictors: all others that have no missing in train
      preds <- setdiff(names(X_train_comp), j)
      dat <- cbind(y = X_train_comp[[j]], X_train_comp[, preds, drop = FALSE])
      dat <- dat[complete.cases(dat), , drop = FALSE]
      if (nrow(dat) < 10) {
        # fallback to mean
        mu <- mean(dat$y, na.rm = TRUE)
        X_test_imp[nas, j] <- mu
      } else {
        # small random forest regressor
        rf_mod <- try(
          randomForest::randomForest(y ~ ., data = dat, ntree = 100),
          silent = TRUE
        )
        if (inherits(rf_mod, "try-error")) {
          mu <- mean(dat$y, na.rm = TRUE)
          X_test_imp[nas, j] <- mu
        } else {
          preds_test <- X_test_imp[nas, preds, drop = FALSE]
          preds_test <- preds_test[, colnames(dat)[-1], drop = FALSE]
          preds_test[is.na(preds_test)] <- 0
          X_test_imp[nas, j] <- predict(rf_mod, preds_test)
        }
      }
    }
  }
  
  as.data.frame(X_test_imp)
}
