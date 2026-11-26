# Changelog

All notable changes to `dimvR` are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.1.2] - 2026-09-16

### Added

- **Downstream model backends**: `evaluate_downstream()` now supports `ranger`, `randomForest`, `glmnet` (plus existing `lm`, `glm`, `xgboost`, `custom`).
- **Runtime dependency checks**: All optional packages validated with actionable errors (`xgboost`, `fastshap`, `future`, `future.apply`, `mice`, `ranger`, `randomForest`, `glmnet`).
- **Test coverage**: 72.83% across 14 test files, 57 test cases (verified 2026-09-14, re-verified 2026-09-16).
- **Test suites added**:
  - `test-dimv-impute-multiple.R`: Multiple imputation structure, reproducibility, edge cases
  - `test-evaluate-downstream-backends.R`: Backend compatibility, missing-package errors, custom models
  - `test-regression-pipeline-edge-cases.R`: Extreme missing rates, multiple mechanisms, matrix input
  - `test-optional-dependency-runtime-checks.R`: Mocked missing-package errors
  - `test-utility-branches.R`: Pooling, missingness, mean imputer, feature-selection fallbacks
- **Progress tracking**: `PROGRESS.md` as single source of truth; `eval/ci_progress_metrics.csv` and `eval/ci_progress_summary.md` regenerated.
- **Coverage gate**: Interim threshold at 40% (`MIN_COVERAGE` env var) enforced in CI.

### Changed

- **`run_full_pipeline()`**: Updated for current `xgboost` interface; re-tested.
- **`dimv_impute_multiple()`**: Validates `m` as positive integer; error message improved.
- **`simulate_missing()`**: Explicit error for unknown mechanisms.
- **Mean imputer**: Warns on non-numeric columns instead of silent failure.
- **Parallel SHAP** (`compute_shap_parallel()`): Graceful fallback to sequential on parallel failure.
- **Documentation**: `README.md` progress snapshot updated to 2026-09-14.

### Fixed

- Multiple imputation preserves observed values exactly (only adds noise to originally missing entries).
- Feature selection handles single-column edge case (returns empty selection).
- MICE backend `apply` validates `k` index bounds.
- Cache helpers (`cache_key`, `get_cached`, `set_cached`) work consistently.

### Experimental (API may change)

- **Feature selection**: `select_features_adaptive()` with methods `adaptive`, `fixed`, `mi`, `hybrid`. Exported but `@keywords internal` in roxygen; output contract not yet stable.
- **Explainability pipeline**: `run_full_pipeline()`, `compute_shap_parallel()`, `generate_report()`.
- **MICE backend**: `mice_backend_train()`, `mice_backend_apply()` — internal, not exported.

---

## [0.1.1] - 2026-04-01

### Added

- **Package rename**: `dimvExplainRpaper` → `dimvR` (DESCRIPTION, NAMESPACE, Rproj).
- **Dependency restructuring**:
  - `Imports`: `MASS`, `stats`, `digest` (minimal core)
  - `Suggests`: `dplyr`, `fastshap`, `future`, `future.apply`, `ggplot2`, `knitr`, `mice`, `ranger`, `randomForest`, `rmarkdown`, `viridis`, `xgboost`, `glmnet`
- **Core imputation API stabilized**:
  - `dimv_train()` — ridge-conditional training with optional feature selection & adaptive lambda
  - `dimv_impute_new()` — deterministic single imputation
  - `dimv_impute_multiple()` — multiple imputation for Rubin pooling
  - `dimv_diagnostics()` — convergence, residual variance, missing RMSE
  - `dimv_convergence_diag()` — iteration, lambda, residual variance summary
  - `adaptive_lambda()` — heuristic 3-level lambda recommendation (exported)
- **Experimental feature selection** (Phase 1):
  - `select_features_adaptive()` — 4 methods: adaptive, fixed, mi, hybrid
  - Internal helpers: `select_adaptive_threshold()`, `select_fixed_threshold()`, `select_mutual_information()`, `select_hybrid()`, `compute_simple_mi()`
- **Experimental explainability pipeline**:
  - `run_full_pipeline()` — end-to-end MCAR/MAR/MNAR benchmark with Rubin pooling
  - `compute_shap_parallel()` — parallel SHAP with `fastshap` + `future`
  - `generate_report()` — HTML report with pooled SHAP MSE, 95% CIs
- **MICE backend abstraction** (internal):
  - `mice_backend_train()` — train-only MICE fitting
  - `mice_backend_apply()` — apply to test set without leakage
- **Mean imputer baseline** (internal): `fit_mean_imputer()`, `apply_mean_imputer()`
- **Diagnostics utilities**: `feature_select_score()`, `plot_feature_selection()`, `dimv_convergence_diag()`, `print_dimv_diag()`
- **Test harness**: `testthat` edition 3, 14 test files
- **CI**: GitHub Actions matrix (Ubuntu/Windows × R release)
- **Smoke benchmark**: Automated iris benchmark in CI
- **Coverage tracking**: `covr` integration, `eval/` artifacts

### Changed

- **Documentation overhaul**: `README.md`, vignette, examples aligned with implementation.
- **Experimental status** explicitly marked in code, docs, and `PROGRESS.md`.

### Fixed

- Namespace imports cleaned (`S3method(print, dimv_imputer)`, explicit `importFrom`).
- Vignette builds without errors.

---

## [0.1.0] - 2025-12-01

### Added

- **Initial DIMV implementation** (R-native):
  - `dimv_train()` — iterative ridge-conditional expectation imputation
  - `dimv_impute_new()` — deterministic imputation on new data
  - `dimv_impute_multiple()` — multiple imputation with Gaussian residual noise
- **Convergence diagnostics**: `diff_history`, residual variances, iteration count
- **Adaptive regularization** (internal): `adaptive_lambda_fn()` based on condition number
- **Feature selection fallback** (internal): correlation-based selector
- **Basic print method**: `print.dimv_imputer()`
- **Reference implementation** of Nguyen et al. (2023) [arXiv:2302.00911]

### Infrastructure

- R package structure with `roxygen2`, `testthat`, `devtools`
- MIT License

---

## Upcoming (Planned for 0.1.3 / 0.2.0)

### Targeted for Q3 2026

- [ ] Raise test coverage from 72.83% → 80% (prioritize XGBoost/SHAP paths)
- [ ] Stabilize `select_features_adaptive()` output contract (field names, types, fallback behavior)
- [ ] Integrate `evaluate_downstream()` into `run_full_pipeline()` for model-agnostic backend support
- [ ] Add non-regression task support (classification) to `run_full_pipeline()`
- [ ] Assess LightGBM/CatBoost for optional backend set
- [ ] Decide on MICE backend: export & stabilize or remove
- [ ] Finalize SHAP/report-generation experimental interfaces
- [ ] CRAN readiness: formal `NEWS.md` (this file), release notes, documentation polish, stricter checks

### Under Consideration

- Kernel-based MI estimator to replace histogram-based `compute_simple_mi()`
- Native parallelization for `dimv_train()` (currently sequential per variable)
- S3 generics for imputer objects (`predict()`, `fitted()`)
- Cross-validation helper for lambda selection

---

## Links

- [GitHub Repository](https://github.com/yht/dimvR)
- [Issues](https://github.com/yht/dimvR/issues)
- [Reference: Nguyen et al. (2023)](https://arxiv.org/abs/2302.00911)
- [Reference: Nguyen et al. (2024)](https://arxiv.org/abs/2407.00411)