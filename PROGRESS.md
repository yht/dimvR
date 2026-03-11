# dimvR Progress Tracker

This document replaces `PROGRESS_PHASE1.md` and summarizes the implementation status that remains relevant to the current state of the repository.

## Current Snapshot

- Package: `dimvR`
- Current focus: core imputation stabilization, documentation alignment, and clearer separation of experimental components
- Latest CI artifact snapshot (2026-02-09):
  - Test files: 8
  - Test cases: 25
  - Estimated coverage: 59.32%
  - Coverage gate: pass at the interim 40% threshold

## Roadmap Execution Notes

### This Week

- [x] Align `README.md` with the current implementation
- [x] Clarify experimental status for:
  - [x] feature selection
  - [x] explainability and benchmarking pipeline
  - [x] internal MICE backend
- [x] Perform an initial `Imports` versus `Suggests` audit
- [x] Reduce the required dependency set in `DESCRIPTION`
- [ ] Remove any remaining suggested dependencies that are no longer directly used
- [ ] Reconcile progress and coverage artifacts so reported metrics are consistent

### This Month

- [ ] Add tests for `dimv_impute_multiple()`
- [ ] Add tests for optional dependency fallbacks
- [ ] Add regression pipeline edge-case tests
- [ ] Reduce coupling between core imputation and explainability/reporting modules
- [ ] Improve installation ergonomics for users who only need the core imputation workflow

### This Quarter

- [ ] Generalize `run_full_pipeline()` beyond regression
- [ ] Stabilize experimental APIs that are candidates for long-term public support
- [ ] Evaluate alternative model backends beyond `xgboost`
- [ ] Improve production and CRAN readiness

## Function and Component Status

### Core Imputation

- [x] `dimv_train()`
- [x] `dimv_impute_new()`
- [x] `dimv_impute_multiple()`
- [x] `dimv_diagnostics()`
- [x] `dimv_convergence_diag()`
- [x] `adaptive_lambda()`

### Feature Selection

- [x] `select_features_adaptive()`
- [x] `select_adaptive_threshold()`
- [x] `select_fixed_threshold()`
- [x] `select_mutual_information()`
- [x] `select_hybrid()`
- [x] `compute_simple_mi()`
- [x] `plot_feature_selection()`
- [ ] Stabilize the feature selection API
- [ ] Expand best-practice documentation

### Explainability and Benchmarking

- [x] `run_full_pipeline()`
- [x] `compute_shap_parallel()`
- [x] `generate_report()`
- [x] CI smoke benchmark workflow
- [ ] Add support for non-`xgboost` backends
- [ ] Add support for non-regression tasks
- [ ] Strengthen optional dependency fallback behavior

### Additional Backends

- [x] Mean imputer helpers
- [x] Internal MICE backend available
- [ ] Export and stabilize the MICE backend
- [ ] Validate benchmarks on additional datasets

## Documentation Checklist

- [x] `README.md` updated to match the current codebase
- [x] Feature selection vignette available
- [x] Example scripts available in `examples/`
- [ ] Add public-facing documentation for experimental workflow limitations
- [ ] Add a formal changelog or release notes

## Infrastructure Checklist

- [x] `testthat` configured and active
- [x] CI checks for Ubuntu and Windows
- [x] Smoke benchmark workflow
- [x] Interim coverage gate
- [ ] Consolidate a single source of truth for coverage metrics
- [ ] Review dependency metadata after the current refactor

## Current Working Assumptions

- Core imputation is treated as the most stable part of the package.
- Feature selection, report generation, SHAP benchmarking, and the internal MICE backend are treated as experimental components.
- Short-term work is prioritized toward documentation alignment, metadata cleanup, and test quality rather than major architectural refactoring.

## Next Actions

1. Remove suggested dependencies that are no longer directly required.
2. Add tests for multiple imputation and optional dependency fallbacks.
3. Reconcile coverage and progress reporting artifacts in CI.
4. Decide which experimental components should be stabilized first.
