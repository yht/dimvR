# dimvR News

## dimvR 0.1.2 (2026-09-16)

### Major Changes

* **Coverage milestone**: Achieved 72.83% test coverage across 14 test files and 57 test cases (verified 2026-09-14, re-verified 2026-09-16). Interim coverage gate (40%) passes in CI.

* **Downstream model backends expanded**: `evaluate_downstream()` now supports `ranger`, `randomForest`, and `glmnet` in addition to `lm`, `glm`, `xgboost`, and custom pre-trained models.

* **Optional dependency fallbacks hardened**: Runtime checks for `xgboost`, `fastshap`, `future`, `future.apply`, `mice`, `ranger`, `randomForest`, `glmnet` with actionable error messages when packages are missing.

* **Test suite expanded**: Added tests for `dimv_impute_multiple()`, regression pipeline edge cases, optional dependency fallbacks, and utility branch coverage.

### New Features

* **Feature selection module** (experimental): `select_features_adaptive()` with four methods — `adaptive` (median-based correlation threshold), `fixed` (user-specified correlation threshold), `mi` (histogram-based mutual information), and `hybrid` (combined correlation + MI). Exported but marked experimental; output contract not yet stable.

* **Diagnostics enhancements**: `feature_select_score()`, `plot_feature_selection()`, `adaptive_lambda()` (heuristic 3-level recommendation), `dimv_convergence_diag()`, `print_dimv_diag()` for model health checks.

* **Multiple imputation**: `dimv_impute_multiple()` generates `m` multiply-imputed datasets with Gaussian residual noise for Rubin-style variance pooling.

* **Report generation**: `generate_report()` creates HTML reports with pooled SHAP MSE, 95% CIs, and summary tables from `run_full_pipeline()` results.

* **MICE backend abstraction** (internal, experimental): `mice_backend_train()` and `mice_backend_apply()` for train-only MICE fitting and test-set application without leakage.

### Infrastructure

* **CI matrix**: Ubuntu and Windows with R release (`.github/workflows/r-check.yml`).

* **Smoke benchmark workflow**: Automated benchmark on iris dataset with MCAR 20% missingness (`.github/workflows/smoke-benchmark.yml`).

* **Coverage gate**: Interim threshold at 40% (`eval/check_coverage_gate.R` with `MIN_COVERAGE` env var).

* **Progress tracking**: `PROGRESS.md` as single source of truth for roadmap and component status.

### Bug Fixes

* `run_full_pipeline()` updated for current `xgboost` interface and re-tested.

* `dimv_impute_multiple()` validates `m` as positive integer.

* `simulate_missing()` rejects unknown mechanisms with clear error.

* Mean imputer helpers handle non-numeric columns with warning.

* Parallel SHAP computation falls back to sequential on failure.

### Deprecated / Experimental Status Clarified

* Feature selection (`select_features_adaptive()` and helpers): explicitly experimental, API may change.

* Explainability & benchmarking pipeline (`run_full_pipeline()`, `compute_shap_parallel()`, `generate_report()`): explicitly experimental.

* Internal MICE backend (`mice_backend_train()`, `mice_backend_apply()`): experimental, not exported.

---

## dimvR 0.1.1 (2026-04-01)

### Major Changes

* **Package renamed** from `dimvExplainRpaper` to `dimvR`.

* **Dependency audit**: Reduced required `Imports` to `MASS`, `stats`, `digest`. Moved `dplyr`, `fastshap`, `future`, `future.apply`, `ggplot2`, `knitr`, `mice`, `ranger`, `randomForest`, `rmarkdown`, `viridis`, `xgboost`, `glmnet` to `Suggests`.

* **Core imputation stabilized**: `dimv_train()`, `dimv_impute_new()`, `dimv_impute_multiple()`, `dimv_diagnostics()`, `dimv_convergence_diag()`, `adaptive_lambda()` marked as stable core.

* **Documentation aligned**: `README.md`, vignettes, and example scripts updated to match implementation boundaries.

* **Test harness**: `testthat` edition 3 configured with 14 test files covering core imputation, diagnostics, multiple imputation, feature selection, experiment pipeline, and optional dependencies.

### New Features

* Initial experimental feature selection implementation (4 methods).

* Initial experimental explainability pipeline with SHAP and Rubin pooling.

* Mean imputer baseline (`fit_mean_imputer()`, `apply_mean_imputer()`).

### Infrastructure

* GitHub Actions CI with R CMD check.

* Automated smoke benchmark in CI.

* Coverage tracking with `covr`.

---

## dimvR 0.1.0 (2025-12-01)

### Initial Release

* **Core DIMV implementation**: R-native conditional ridge regression imputation with iterative convergence (`dimv_train()`, `dimv_impute_new()`).

* **Deterministic and multiple imputation**: Single imputation and `m`-fold multiple imputation with Gaussian residual noise.

* **Convergence diagnostics**: Iteration tracking, diff history, residual variance summaries.

* **Adaptive regularization**: Per-variable lambda based on condition number (internal).

* **Basic feature selection**: Correlation-based fallback selector (internal).

* **Reference**: Nguyen et al. (2023) "Conditional expectation with regularization for missing data imputation" (arXiv:2302.00911).