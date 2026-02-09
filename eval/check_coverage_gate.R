#!/usr/bin/env Rscript

metrics_path <- file.path("eval", "ci_progress_metrics.csv")
if (!file.exists(metrics_path)) {
  stop("Coverage gate: metrics file not found: ", metrics_path)
}

metrics <- read.csv(metrics_path, stringsAsFactors = FALSE)
if (!all(c("metric", "value") %in% names(metrics))) {
  stop("Coverage gate: invalid metrics file format.")
}

idx <- which(metrics$metric == "coverage_percent")
if (length(idx) != 1) {
  stop("Coverage gate: coverage_percent metric not found.")
}

coverage_val <- suppressWarnings(as.numeric(metrics$value[idx]))
if (!is.finite(coverage_val)) {
  stop("Coverage gate failed: coverage_percent is not numeric (NA/NaN/Inf).")
}

min_cov <- Sys.getenv("MIN_COVERAGE", unset = "60")
min_cov_num <- suppressWarnings(as.numeric(min_cov))
if (!is.finite(min_cov_num)) {
  stop("Coverage gate: MIN_COVERAGE is not numeric: ", min_cov)
}

cat(sprintf("Coverage gate: coverage_percent=%.2f, threshold=%.2f\n", coverage_val, min_cov_num))

summary_file <- file.path("eval", "ci_coverage_gate_summary.md")
summary_lines <- c(
  "# Coverage Gate",
  "",
  sprintf("- coverage_percent: %.2f", coverage_val),
  sprintf("- threshold (MIN_COVERAGE): %.2f", min_cov_num),
  "",
  if (coverage_val >= min_cov_num) "Result: PASS" else "Result: FAIL"
)
writeLines(summary_lines, con = summary_file)

if (coverage_val < min_cov_num) {
  stop(sprintf("Coverage gate failed: %.2f < %.2f", coverage_val, min_cov_num))
}
