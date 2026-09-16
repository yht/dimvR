# dimvR Progress Tracker

This document replaces `PROGRESS_PHASE1.md` and summarizes the implementation status that remains relevant to the current state of the repository.

## Current Snapshot (Q3 2026)
- Package: `dimvR` (v0.1.2)
- Current focus: core imputation stabilization, documentation alignment, and clearer separation of experimental components
- Coverage status Q3 2026: **72.83%** (measured 2026-09-14 with `covr`; latest checked-in artifact)
- Q3 2026 status: Core imputation is functional; experimental features remain under active refinement
- Latest coverage verification (2026-09-14), with the test suite re-run on 2026-09-16:
  - Test files: 14 (in `tests/testthat/`)
  - Test cases: 57
  - Coverage percent: **72.83%**
  - Coverage gate: pass at interim threshold (40%); Q3 target 80% not yet verified
  - Coverage run completed with 8 optional-dependency tests skipped
  - Core imputation: fully functional (100% of core items)
  - Feature selection: 7/9 items complete (77.8%), still EXPERIMENTAL
  - MICE backend: 2/4 items complete (50%), experimental, not exported
  - Downstream evaluator: supports `lm`, `glm`, `ranger`, `randomForest`, `glmnet`, `xgboost`, and custom models
  - End-to-end xgboost pipeline: functional for regression workflows only

### Q3 2026 Priority Tags
- `🟢` = Complete & stable
- `🟡` = Experimental, usable with cautions
- `🔴` = Not yet implemented / requires further work

## Roadmap Execution Notes

### Completed Through 2026-04-01

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
- [x] Add tests for optional dependency fallbacks

### Remaining For Q3 2026

- [ ] Raise automated test coverage from 72.83% toward the 80% target, prioritizing skipped XGBoost/SHAP paths
- [ ] Finalize the stable contract for `select_features_adaptive()` and expand edge-case tests
- [ ] Integrate `evaluate_downstream()` into `run_full_pipeline()` and add non-regression task support
- [x] Evaluate alternative downstream model backends beyond `xgboost` (`ranger`, `randomForest`, and `glmnet`)
- [ ] Assess LightGBM/CatBoost and decide whether they belong in the supported optional-backend set
- [ ] Decide whether the MICE backend should become a stable exported API
- [ ] Finalize SHAP and report-generation experimental interfaces
- [ ] Complete CRAN-readiness checks, dependency review, documentation polish, and release notes

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

- Completion: 85.7% (6 of 7 items completed)
- [x] `run_full_pipeline()`
- [x] `compute_shap_parallel()`
- [x] `generate_report()`
- [x] CI smoke benchmark workflow
- [x] Strengthen optional dependency fallback behavior
- [x] Add non-`xgboost` backends to `evaluate_downstream()`
- [ ] Add support for non-regression tasks

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

- Completion: 85.7% (6 of 7 items completed)
- [x] `testthat` configured and active
- [x] CI checks for Ubuntu and Windows
- [x] Smoke benchmark workflow
- [x] Interim coverage gate
- [x] Regenerated local progress and coverage artifacts against the current repository state (72.83% coverage)
- [x] Treat the root `eval/` artifacts as the active coverage/progress baseline
- [x] Consolidate a single source of truth for coverage metrics
- [ ] Review dependency metadata after the current refactor

## Current Working Assumptions

- Core imputation is treated as the most stable part of the package.
- Feature selection, report generation, SHAP benchmarking, and the internal MICE backend are treated as experimental components.
- Short-term work is prioritized toward documentation alignment, metadata cleanup, and test quality rather than major architectural refactoring.

## Priorities for 2026-09-16

1. Raise coverage from 72.83% toward 80%, especially for skipped XGBoost/SHAP paths.
2. Finalize and document the `select_features_adaptive()` output contract.
3. Connect the model-agnostic `evaluate_downstream()` helper to the end-to-end experiment workflow.
4. Review optional dependency metadata and make the experimental API boundaries explicit.
5. Complete CRAN-readiness checks and add formal release notes.
