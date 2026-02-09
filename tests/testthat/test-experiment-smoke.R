test_that("simulate_missing supports MCAR/MAR/MNAR and rejects unknown mechanism", {
  set.seed(101)
  X <- as.data.frame(matrix(rnorm(120), ncol = 4))

  for (mech in c("MCAR", "MAR", "MNAR")) {
    Xm <- dimvR:::simulate_missing(X, rate = 0.2, mechanism = mech)
    expect_equal(dim(Xm), dim(X))
    expect_true(anyNA(Xm))
  }

  expect_error(
    dimvR:::simulate_missing(X, rate = 0.2, mechanism = "UNKNOWN"),
    "Unknown mechanism"
  )
})

test_that("run_full_pipeline smoke test returns expected columns", {
  skip_if_not_installed("xgboost")
  skip_if_not_installed("fastshap")
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  set.seed(123)
  X <- data.frame(
    x1 = rnorm(36),
    x2 = rnorm(36),
    x3 = rnorm(36)
  )
  y <- 1.2 * X$x1 - 0.5 * X$x2 + rnorm(36, sd = 0.2)

  # Guard against interface differences in newer xgboost releases.
  probe_ok <- tryCatch({
    mod_probe <- suppressWarnings(xgboost::xgboost(
      x = as.matrix(X),
      y = y,
      nrounds = 2,
      objective = "reg:squarederror",
      verbose = 0
    ))
    pred_probe <- predict(mod_probe, xgboost::xgb.DMatrix(as.matrix(X[1:2, , drop = FALSE])))
    is.numeric(pred_probe)
  }, error = function(e) FALSE)
  if (!isTRUE(probe_ok)) {
    skip("xgboost interface incompatible with current run_full_pipeline implementation")
  }

  out <- suppressWarnings(run_full_pipeline(
    X = X,
    y = y,
    missing_rates = 0.1,
    mechanisms = "MCAR",
    imputers = "mean",
    nsim = 5,
    workers = 1,
    m_imp = 2,
    seed = 2026
  ))

  expect_true(is.data.frame(out))
  expect_true(nrow(out) >= 1)
  expect_true(all(c("method", "mechanism", "rate", "mse_pred", "mse_shap", "m") %in% names(out)))
})
