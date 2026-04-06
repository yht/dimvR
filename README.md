# dimvR

> Regularized Conditional Distribution-based Imputation for Missing Data in R

[![R-CMD-check](https://github.com/yht/dimvR/actions/workflows/r-check.yml/badge.svg)](https://github.com/yht/dimvR/actions/workflows/r-check.yml)
[![Smoke Benchmark](https://github.com/yht/dimvR/actions/workflows/smoke-benchmark.yml/badge.svg)](https://github.com/yht/dimvR/actions/workflows/smoke-benchmark.yml)

This package implements a regularized conditional distribution-based approach to missing data imputation,
commonly referred to as **DIMV (Distribution-based Imputation using Conditional Expectation with Regularization)**
[(Nguyen et al., 2023)](https://arxiv.org/abs/2302.00911). The method imputes missing entries by estimating the
conditional expectation of each variable given the others while applying regularization to improve numerical
stability and robustness in high-dimensional or multicollinear settings.

The package is motivated by methodological advances discussed in
*“Explainability of Machine Learning Models under Missing Data”* [(Nguyen et al., 2024)](https://arxiv.org/abs/2407.00411v1),
which emphasizes the importance of principled imputation for preserving both predictive performance and model
interpretability. The **dimvR** package provides a reproducible R implementation of the
[DIMV framework](https://github.com/maianhpuco/DIMVImputation) and currently supports core imputation workflows,
with additional experimental components for benchmarking, explainability, and reporting.

## Current Capability Summary

- **Regularized Conditional Imputation (DIMV):**  
  Implements conditional expectation-based imputation with ridge regularization for enhanced stability in
  high-dimensional or multicollinear settings.

- **Deterministic and Multiple Imputation:**  
  Supports both deterministic imputation and multiple imputation through Gaussian residual noise injection to
  approximate uncertainty propagation.

- **Convergence-Based Iterative Procedure:**  
  Iteratively updates missing values variable-by-variable until convergence to support stable and internally
  consistent imputations.

- **Feature Selection (Experimental):**  
  Provides adaptive, fixed-threshold, mutual information, and hybrid selection strategies to reduce noise and
  multicollinearity in conditional models. This component remains explicitly experimental.

- **Explainability & Benchmark Pipeline (Experimental):**  
  Provides a regression-oriented evaluation workflow with SHAP-based explainability, pooled benchmarking, and
  HTML report generation. This component remains experimental.

- **R-Native Core Imputation:**  
  Core DIMV training and imputation are implemented in R with a smaller required dependency set than the
  experimental benchmarking and reporting workflows.

- **Downstream Model Compatibility:**  
  Imputed outputs can be used in downstream statistical or machine learning workflows. The bundled experimental
  pipeline currently targets regression using `xgboost`.

## Engineering Progress (P0)

| Item | Status | Evidence |
|---|---|---|
| CI matrix (Ubuntu/Windows, R release) | Implemented | `.github/workflows/r-check.yml` |
| Standard test harness (`testthat`) | Implemented | `tests/testthat.R`, `DESCRIPTION` (`Config/testthat/edition: 3`) |
| Automated smoke benchmark in CI | Implemented | `.github/workflows/smoke-benchmark.yml`, `eval/smoke_benchmark.R` |
| Benchmark artifact publication | Implemented | CI artifact: `ci_smoke_benchmark.csv`, `ci_smoke_summary.md` |
| Progress metrics (test count + coverage) | Implemented | `eval/ci_progress_metrics.csv`, `eval/ci_progress_summary.md` |
| Interim coverage gate | Implemented | `MIN_COVERAGE=40` via `eval/check_coverage_gate.R` in CI workflow |

Execution snapshot:
- Verified active snapshot (2026-04-06): 11 test files, 34 test cases, 76.87% estimated coverage
- Coverage gate status at that verified snapshot: pass against interim threshold (40%)
- Tracked implementation completion: 70.4% (based on `PROGRESS.md`)

Current implementation notes:
- `run_full_pipeline()` is currently a regression-oriented workflow using `xgboost`.
- Feature selection helpers are exported but explicitly marked experimental in `R/feature_selection.R`.
- The internal MICE backend exists in `R/mice_backend.R` but is still marked experimental and is not exported.
- Core imputation uses a smaller required dependency set; benchmarking, SHAP, and report-generation helpers rely on additional suggested packages.
- `dimvExplainR/eval/` is the active checked-in progress/coverage baseline; the root `eval/` folder is archival unless explicitly regenerated.

## Feature Selection Direction

Current audit conclusion:
- `select_features_adaptive()` is already the de facto public entry point for feature selection because it is exported, tested, documented, and used by `dimv_train()`.
- The current output shape is useful but not yet a stable API contract; field meanings remain partly ambiguous across `adaptive`, `fixed`, `mi`, and `hybrid` modes.
- Internal helper functions such as `select_adaptive_threshold()`, `select_fixed_threshold()`, `select_mutual_information()`, `select_hybrid()`, and `compute_simple_mi()` should remain implementation details rather than becoming user-facing stability commitments.

Proposed stable contract for `select_features_adaptive()`:
- Keep `selected_features` as the primary downstream field and guarantee that it is always an integer vector of valid column indexes.
- Keep `selected_names` aligned 1:1 with `selected_features`.
- Normalize target metadata into stable fields such as `target_index` and `target_name`.
- Add explicit metadata for `method`, `n_selected`, `candidate_count`, `threshold_used`, and `fallback_used`.
- Replace the current ambiguous score semantics with a stable ranking object that can carry correlation, MI, hybrid, and final scores without changing the top-level contract.

Immediate stabilization to-do:
- Tighten argument validation for `min_features`, `max_features`, `threshold`, `nbins`, and invalid `target_var` cases.
- Make fallback behavior explicit in the returned object instead of silently switching to top-correlation selection.
- Add contract-oriented tests that validate output structure and edge-case behavior, not just successful execution.
- Update vignette and reference docs so the experimental label remains, but the intended stable interface boundaries are documented clearly.


## Current Scope, Assumptions, and Limitations

**Assumptions**:

- **Approximate Multivariate Normality:**  
  DIMV implicitly relies on the assumption that the joint distribution of variables is approximately multivariate
  normal, enabling the use of conditional expectation as the optimal estimator.

- **Missing at Random (MAR) Conditions Preferred:**  
  The method performs best when data are Missing Completely at Random (MCAR) or Missing at Random (MAR).  
  Under Missing Not at Random (MNAR), bias may persist unless additional modeling is applied.

- **Linear Conditional Relationships:**  
  Ridge-based conditional models assume linear relationships among variables. Although robust under moderate
  deviations, highly nonlinear dependencies may reduce accuracy.

**Limitations**:

- **Not Tailored for Strong Nonlinear Interactions:**  
  In datasets where nonlinearities dominate, tree-based or neural-network-based imputation may outperform DIMV.

- **Performance Degradation in High Missingness (>60%):**  
  When a large proportion of entries are missing, conditional estimates become less reliable, especially for 
  variables with low correlation to others.

- **Multiple Imputation is Approximate:**  
  Noise injection relies on Gaussian residual estimates, which may inadequately capture uncertainty if residuals 
  deviate significantly from normality.

- **Experimental Pipeline Components:**  
  Feature selection, SHAP benchmarking, report generation, and the internal MICE backend are still under active
  refinement and should be treated as experimental interfaces.

- **Current Experiment Scope is Regression-Focused:**  
  The bundled end-to-end experiment runner currently targets regression with `xgboost`, so classification and
  alternative model backends are not yet first-class workflow options.


## References

Vu, M. A., Nguyen, T., Do, T. T., Phan, N., Chawla, N. V., Halvorsen, P., Riegler, M. A. & Nguyen, B. T. (2023).
**Conditional expectation with regularization for missing data imputation**. *arXiv:2302.00911v3 [stat.ML]*.  
https://arxiv.org/abs/2302.00911

Nguyen, M.-H., Dao, M.-H., Tran, M.-T., & Phung, D. (2024).  
**Explainability of Machine Learning Models under Missing Data**. *arXiv:2407.00411v3 [cs.LG]*.  
https://arxiv.org/abs/2407.00411

Van Buuren, S., & Groothuis-Oudshoorn, K. (2011).  
**mice: Multivariate Imputation by Chained Equations in R**. *Journal of Statistical Software*, 45(3).  
https://doi.org/10.18637/jss.v045.i03

Rubin, D. B. (1987).  
**Multiple Imputation for Nonresponse in Surveys**. *John Wiley & Sons*.
https://onlinelibrary.wiley.com/doi/book/10.1002/9780470316696

Little, R. J. A., & Rubin, D. B. (2019).  
**Statistical Analysis with Missing Data (3rd ed.)**. *Wiley*.
https://onlinelibrary.wiley.com/doi/book/10.1002/9781119482260


## Example Usage

```r
library(dimvR)

set.seed(123)

# Simulated dataset with missing values
X <- data.frame(
  x1 = rnorm(100),
  x2 = rnorm(100),
  x3 = rnorm(100)
)

# Introduce missingness
for (j in seq_along(X)) {
  X[sample(1:nrow(X), 10), j] <- NA
}

# Train DIMV imputer
imp <- dimv_train(X, lambda = 0.1, feature_select = TRUE)

# Print summary
print(imp)

# Impute new data
X_imputed <- dimv_impute_new(imp, X)
head(X_imputed)

# Multiple imputation (m = 5)
imputations <- dimv_impute_multiple(imp, X, m = 5)

# Feature selection (standalone)
fs <- select_features_adaptive(X, target_var = "x1", method = "hybrid", max_features = 2)
fs$selected_features
```

## Documentation and Examples

- Vignette: `vignettes/feature_selection.Rmd`
- Examples: `examples/day1_feature_selection_demo.R`, `examples/benchmark_feature_selection.R`


## Roadmap Overview

### Completed Through 2026-04-01

- Align package metadata and documentation with the current implementation boundaries.
- Clarify experimental status for feature selection, SHAP benchmarking, report generation, and the internal MICE backend.
- Audit and simplify `Imports` versus `Suggests` so optional workflows do not overstate core install requirements.
- Remove stale optional metadata/runtime checks that no longer match the active implementation.
- Add direct tests for `dimv_impute_multiple()`.
- Add regression pipeline edge-case tests.
- Add runtime-check coverage tests for optional dependency paths.

### Planned For April 2026

- Add tests for optional dependency fallback behavior.
- Stabilize the `select_features_adaptive()` contract before expanding feature-selection scope.
- Reduce coupling between core imputation and explainability/reporting modules.
- Improve installation and dependency ergonomics for users who only need the core imputation workflow.
- Build on the refreshed 2026-04-06 coverage snapshot with additional low-risk tests to move toward the April 78% checkpoint.

### Planned For Q2 2026

- Generalize the experiment pipeline beyond regression-only evaluation.
- Stabilize experimental APIs and decide which components should become long-term public interfaces.
- Evaluate alternative model backends beyond `xgboost`.
- Advance production and CRAN readiness with tighter checks, documentation polish, and stronger quality gates.

For more detailed execution tracking, see `PROGRESS.md`.
