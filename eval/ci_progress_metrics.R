#!/usr/bin/env Rscript

test_files <- list.files("tests/testthat", pattern = "^test.*\\.R$", full.names = TRUE)
test_count <- 0L
if (length(test_files) > 0) {
  for (f in test_files) {
    lines <- readLines(f, warn = FALSE)
    test_count <- test_count + sum(grepl("test_that\\s*\\(", lines))
  }
}

coverage_pct <- NA_real_
if (requireNamespace("covr", quietly = TRUE)) {
  cov <- tryCatch(
    covr::package_coverage(type = "tests", quiet = TRUE),
    error = function(e) NULL
  )
  if (!is.null(cov)) {
    coverage_pct <- round(covr::percent_coverage(cov), 2)
  }
}

metrics <- data.frame(
  metric = c("test_files", "test_cases", "coverage_percent"),
  value = c(length(test_files), test_count, coverage_pct),
  stringsAsFactors = FALSE
)

if (!dir.exists("eval")) {
  dir.create("eval", recursive = TRUE)
}

write.csv(metrics, file = file.path("eval", "ci_progress_metrics.csv"), row.names = FALSE)

summary_lines <- c(
  "# CI Progress Metrics",
  "",
  paste0("- Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "| metric | value |",
  "|---|---:|",
  sprintf("| test_files | %d |", length(test_files)),
  sprintf("| test_cases | %d |", test_count),
  sprintf("| coverage_percent | %s |", ifelse(is.na(coverage_pct), "NA", format(coverage_pct, nsmall = 2)))
)

writeLines(summary_lines, con = file.path("eval", "ci_progress_summary.md"))
