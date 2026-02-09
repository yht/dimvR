#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dimvR)
})

set.seed(123)

if (!dir.exists("eval")) {
  dir.create("eval", recursive = TRUE)
}

# Small and fast benchmark intended for CI smoke checks.
data(iris)
X <- iris[, 1:4]
n <- nrow(X)
p <- ncol(X)
mask <- matrix(runif(n * p) < 0.2, nrow = n, ncol = p)
X_miss <- X
X_miss[mask] <- NA

rmse_eval <- function(original, imputed, missing_ref) {
  idx <- is.na(as.matrix(missing_ref))
  diff <- as.matrix(original)[idx] - as.matrix(imputed)[idx]
  sqrt(mean(diff^2))
}

bench_dimv <- system.time({
  imp <- dimv_train(X_miss, maxit = 20, tol = 1e-4)
  X_dimv <- dimv_impute_new(imp, X_miss, maxit = 20, tol = 1e-4)
})

bench_mean <- system.time({
  X_mean <- X_miss
  for (j in seq_len(ncol(X_mean))) {
    v <- X_mean[[j]]
    v[is.na(v)] <- mean(v, na.rm = TRUE)
    X_mean[[j]] <- v
  }
})

results <- data.frame(
  method = c("dimvR", "mean"),
  rmse = c(rmse_eval(X, X_dimv, X_miss), rmse_eval(X, X_mean, X_miss)),
  runtime_sec = c(unname(bench_dimv[["elapsed"]]), unname(bench_mean[["elapsed"]])),
  missing_rate = 0.2,
  dataset = "iris",
  stringsAsFactors = FALSE
)

write.csv(results, file = file.path("eval", "ci_smoke_benchmark.csv"), row.names = FALSE)

summary_lines <- c(
  "# CI Smoke Benchmark",
  "",
  paste0("- Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "- Dataset: iris",
  "- Missingness: MCAR 20%",
  "",
  "## Results",
  "",
  "| method | rmse | runtime_sec |",
  "|---|---:|---:|",
  sprintf("| %s | %.4f | %.4f |", results$method[1], results$rmse[1], results$runtime_sec[1]),
  sprintf("| %s | %.4f | %.4f |", results$method[2], results$rmse[2], results$runtime_sec[2])
)

writeLines(summary_lines, con = file.path("eval", "ci_smoke_summary.md"))
