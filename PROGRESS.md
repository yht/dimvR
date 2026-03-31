# dimvR Progress Tracker

This document replaces `PROGRESS_PHASE1.md` and summarizes the implementation status that remains relevant to the current state of the repository.

## Current Snapshot

- Package: `dimvR`
- Current focus: core imputation stabilization, documentation alignment, and clearer separation of experimental components
- Overall tracked completion: 68.5% (37 of 54 tracked items completed)
- Latest executed local verification snapshot (2026-03-31):
  - Test files: 10
  - Test cases: 29
  - Coverage percent: 76.19
  - Coverage gate: pass at the interim 40% threshold
  - Local `testthat::test_local('.')`: pass

## Roadmap Execution Notes

### Completed Through 2026-03-31

- [x] Align `README.md` with the current implementation
- [x] Clarify experimental status for:
  - [x] feature selection
  - [x] explainability and benchmarking pipeline
  - [x] internal MICE backend
- [x] Perform an initial `Imports` versus `Suggests` audit
- [x] Reduce the required dependency set in `DESCRIPTION`
- [x] Remove any remaining suggested dependencies that are no longer directly used
- [x] Audit progress and coverage artifacts and reconcile the checked-in reporting baseline
- [x] Add tests for `dimv_impute_multiple()`
- [x] Add regression pipeline edge-case tests
- [x] Update `run_full_pipeline()` for the active `xgboost` interface and rerun tests

### Planned For April 2026

- [ ] Add tests for optional dependency fallbacks
- [ ] Reduce coupling between core imputation and explainability/reporting modules
- [ ] Improve installation ergonomics for users who only need the core imputation workflow

### Planned For Q2 2026

- [ ] Generalize `run_full_pipeline()` beyond regression
- [ ] Stabilize experimental APIs that are candidates for long-term public support
- [ ] Evaluate alternative model backends beyond `xgboost`
- [ ] Improve production and CRAN readiness

## Function and Component Status

### Core Imputation

- Completion: 100.0% (6 of 6 items completed)
- [x] `dimv_train()`
- [x] `dimv_impute_new()`
- [x] `dimv_impute_multiple()`
- [x] `dimv_diagnostics()`
- [x] `dimv_convergence_diag()`
- [x] `adaptive_lambda()`

### Feature Selection

- Completion: 77.8% (7 of 9 items completed)
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

- Completion: 57.1% (4 of 7 items completed)
- [x] `run_full_pipeline()`
- [x] `compute_shap_parallel()`
- [x] `generate_report()`
- [x] CI smoke benchmark workflow
- [ ] Add support for non-`xgboost` backends
- [ ] Add support for non-regression tasks
- [ ] Strengthen optional dependency fallback behavior

### Additional Backends

- Completion: 50.0% (2 of 4 items completed)
- [x] Mean imputer helpers
- [x] Internal MICE backend available
- [ ] Export and stabilize the MICE backend
- [ ] Validate benchmarks on additional datasets

## Documentation Checklist

- Completion: 80.0% (4 of 5 items completed)
- [x] `README.md` updated to match the current codebase
- [x] Feature selection vignette available
- [x] Example scripts available in `examples/`
- [x] Add public-facing documentation for experimental workflow limitations
- [ ] Add a formal changelog or release notes

## Infrastructure Checklist

- Completion: 71.4% (5 of 7 items completed)
- [x] `testthat` configured and active
- [x] CI checks for Ubuntu and Windows
- [x] Smoke benchmark workflow
- [x] Interim coverage gate
- [x] Regenerated local progress and coverage artifacts against the current repository state
- [ ] Consolidate a single source of truth for coverage metrics
- [ ] Review dependency metadata after the current refactor

## Current Working Assumptions

- Core imputation is treated as the most stable part of the package.
- Feature selection, report generation, SHAP benchmarking, and the internal MICE backend are treated as experimental components.
- Short-term work is prioritized toward documentation alignment, metadata cleanup, and test quality rather than major architectural refactoring.

## Priorities for 2026-04-01

1. Complete the dependency cleanup review, including any remaining stale `Suggests` entries and runtime checks.
2. Add targeted tests for optional dependency fallback paths.
3. Review CI-generated progress and coverage artifacts and define a single reporting baseline.
4. Identify and document the next small implementation gaps that can be closed without major architectural changes.
