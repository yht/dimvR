test_that(".require_pkg errors for a missing package", {
  expect_error(
    dimvR:::.require_pkg("definitelyNotARealPackage"),
    "required but not installed"
  )
})

test_that("run_full_pipeline requests mice only when the mice imputer is selected", {
  calls <- character()

  testthat::local_mocked_bindings(
    .require_pkg = function(pkg) {
      calls <<- c(calls, pkg)
      if (identical(pkg, "mice")) {
        stop("Package 'mice' is required but not installed.")
      }
      invisible(TRUE)
    },
    .package = "dimvR"
  )

  X <- data.frame(x1 = c(1, 2, 3, 4), x2 = c(4, 3, 2, 1))
  y <- c(1, 2, 3, 4)

  expect_error(
    run_full_pipeline(
      X = X,
      y = y,
      imputers = "mice",
      missing_rates = 0.2,
      mechanisms = "MCAR",
      nsim = 1,
      workers = 1,
      m_imp = 2,
      seed = 99
    ),
    "Package 'mice' is required but not installed."
  )

  expect_equal(calls, c("xgboost", "fastshap", "mice"))
})

test_that("compute_shap_parallel stops when an optional SHAP dependency is unavailable", {
  testthat::local_mocked_bindings(
    requireNamespace = function(pkg, quietly = TRUE) {
      !identical(pkg, "fastshap")
    },
    .package = "base"
  )

  expect_error(
    compute_shap_parallel(
      model = NULL,
      X_train = data.frame(x = 1:3),
      X_test = data.frame(x = 1:2),
      nsim = 1,
      n_workers = 1
    ),
    "Package 'fastshap' is required for compute_shap_parallel\\(\\)"
  )
})

test_that("generate_report stops when a report dependency is unavailable", {
  testthat::local_mocked_bindings(
    requireNamespace = function(pkg, quietly = TRUE) {
      !identical(pkg, "ggplot2")
    },
    .package = "base"
  )

  results <- data.frame(
    method = "mean",
    mechanism = "MCAR",
    rate = 0.2,
    m = 1,
    mse_pred = 0.1,
    mse_pred_se = 0.01,
    mse_shap = 0.2,
    mse_shap_se = 0.02,
    stringsAsFactors = FALSE
  )

  expect_error(
    generate_report(results, out_file = tempfile(fileext = ".html")),
    "Package 'ggplot2' is required but not installed"
  )
})

test_that("plot_feature_selection stops when ggplot2 is unavailable", {
  testthat::local_mocked_bindings(
    requireNamespace = function(pkg, quietly = TRUE) {
      !identical(pkg, "ggplot2")
    },
    .package = "base"
  )

  expect_error(
    plot_feature_selection(data.frame(feature = "x1", score = 1)),
    "Package 'ggplot2' is required for plot_feature_selection\\(\\)"
  )
})
