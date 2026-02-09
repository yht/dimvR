test_that("generate_report creates HTML output from minimal pooled results", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("knitr")
  skip_if_not_installed("viridis")

  results <- data.frame(
    method = c("dimv_R", "mean"),
    mechanism = c("MCAR", "MCAR"),
    rate = c(0.2, 0.2),
    m = c(5, 5),
    mse_pred = c(0.11, 0.25),
    mse_pred_se = c(0.02, 0.03),
    mse_shap = c(0.07, 0.18),
    mse_shap_se = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  out_file <- tempfile(fileext = ".html")
  expect_invisible(generate_report(results, out_file = out_file, dataset_name = "Smoke Dataset"))
  expect_true(file.exists(out_file))

  html <- readLines(out_file, warn = FALSE)
  txt <- paste(html, collapse = "\n")
  expect_match(txt, "DIMV Explainability Report")
  expect_match(txt, "Smoke Dataset")
})
