library(testthat)
library(dimvR)
test_that("run_full_pipeline handles edge cases in missing data simulation", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("fastshap")
  set.seed(42)
  X <- data.frame(
    x1 = rnorm(20),
    x2 = rnorm(20),
    x3 = rnorm(20),
    x4 = rnorm(20),
    x5 = rnorm(20),
    x6 = rnorm(20),
    x7 = rnorm(20),
    x8 = rnorm(20),
    x9 = rnorm(20),
    x10 = rnorm(20)
  )
  y <- rnorm(20)
  expect_silent({
    result <- run_full_pipeline(
      X = X, y = y,
      missing_rates = 0.01,
      mechanisms = "MCAR",
      imputers = "dimv_R",
      nsim = 10,
      workers = 1,
      m_imp = 2,
      seed = 123
    )
    expect_s3_class(result, "data.frame")
    expect_named(result, c("method", "mechanism", "rate", "mse_pred", "mse_pred_se", 
                          "mse_shap", "mse_shap_se", "m"))
    expect_equal(nrow(result), 1)
  })
  expect_silent({
    result <- run_full_pipeline(
      X = X, y = y,
      missing_rates = 0.5,
      mechanisms = "MCAR",
      imputers = "dimv_R",
      nsim = 10,
      workers = 1,
      m_imp = 2,
      seed = 123
    )
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 1)
  })
  expect_silent({
    result <- run_full_pipeline(
      X = X, y = y,
      missing_rates = 0.2,
      mechanisms = c("MCAR", "MAR", "MNAR"),
      imputers = "dimv_R",
      nsim = 10,
      workers = 1,
      m_imp = 2,
      seed = 123
    )
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 3)
  })
})
test_that("run_full_pipeline validates input parameters", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("fastshap")
  set.seed(42)
  X <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  y <- rnorm(20)
  expect_error({
    run_full_pipeline(
      X = X, y = y,
      missing_rates = 0.2,
      mechanisms = "INVALID_MECH",
      imputers = "dimv_R",
      nsim = 10,
      workers = 1,
      m_imp = 2,
      seed = 123
    )
  }, "Unknown mechanism")
  expect_error({
    run_full_pipeline(
      X = X, y = y,
      missing_rates = 0.2,
      mechanisms = "MCAR",
      imputers = "invalid_imputer",
      nsim = 10,
      workers = 1,
      m_imp = 2,
      seed = 123
    )
  }, "Unknown imputer")
})
test_that("run_full_pipeline works with different data types", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("fastshap")
  set.seed(42)
  X_mat <- matrix(rnorm(200), ncol = 10)
  colnames(X_mat) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10")
  y_vec <- rnorm(20)
  expect_silent({
    result <- run_full_pipeline(
      X = X_mat, y = y_vec,
      missing_rates = 0.2,
      mechanisms = "MCAR",
      imputers = "dimv_R",
      nsim = 10,
      workers = 1,
      m_imp = 2,
      seed = 123
    )
    expect_s3_class(result, "data.frame")
    expect_equal(nrow(result), 1)
  })
})
