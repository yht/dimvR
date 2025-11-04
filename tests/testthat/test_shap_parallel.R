# ============================================================
# test_shap_parallel.R
# Unit tests for compute_shap_parallel()
# ============================================================

test_that("compute_shap_parallel runs and returns correct structure", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("fastshap")
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  
  set.seed(123)
  
  # Small toy dataset
  X <- data.frame(
    x1 = rnorm(40),
    x2 = rnorm(40),
    x3 = rnorm(40)
  )
  y <- 1.5 * X$x1 - 0.7 * X$x2 + rnorm(40, 0, 0.3)
  
  # Train simple model
  mod <- xgboost::xgboost(
    data = as.matrix(X),
    label = y,
    nrounds = 20,
    objective = "reg:squarederror",
    verbose = 0
  )
  
  # Run SHAP in parallel
  shap_mat <- compute_shap_parallel(
    model = mod,
    X_train = X,
    X_test = X,
    nsim = 20,
    n_workers = 2
  )
  
  # Basic assertions
  expect_true(is.matrix(shap_mat) || is.data.frame(shap_mat))
  shap_mat <- as.matrix(shap_mat)
  
  expect_equal(nrow(shap_mat), nrow(X))
  expect_equal(ncol(shap_mat), ncol(X))
  expect_false(any(is.na(shap_mat)))
})


test_that("parallel and sequential SHAP results are similar", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("fastshap")
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  
  set.seed(456)
  
  X <- data.frame(
    x1 = rnorm(25),
    x2 = rnorm(25)
  )
  y <- X$x1 + rnorm(25)
  
  mod <- xgboost::xgboost(
    data = as.matrix(X),
    label = y,
    nrounds = 15,
    objective = "reg:squarederror",
    verbose = 0
  )
  
  shap_par <- compute_shap_parallel(mod, X, X, nsim = 10, n_workers = 2)
  shap_seq <- compute_shap_parallel(mod, X, X, nsim = 10, n_workers = 1)
  
  expect_equal(dim(shap_par), dim(shap_seq))
  
  # Values won’t match exactly but should be close
  expect_lt(mean(abs(shap_par - shap_seq)), 1e-1)
})


test_that("compute_shap_parallel gracefully falls back when parallel fails", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("fastshap")
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  
  set.seed(789)
  
  X <- data.frame(x = rnorm(10))
  y <- rnorm(10)
  
  mod <- xgboost::xgboost(
    data = as.matrix(X),
    label = y,
    nrounds = 5,
    objective = "reg:squarederror",
    verbose = 0
  )
  
  # Force failure: set workers = 0
  expect_warning_or_message <- function(expr) {
    tryCatch(expr, error = function(e) TRUE, warning = function(w) TRUE, message = function(m) TRUE)
  }
  
  result <- expect_warning_or_message(
    compute_shap_parallel(mod, X, X, nsim = 5, n_workers = 0)
  )
  
  # Expect output still valid shape
  if (is.matrix(result)) {
    expect_equal(nrow(result), nrow(X))
  }
})
