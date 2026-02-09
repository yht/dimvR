# dimvR

> Regularized Conditional Distribution-based Imputation for Missing Data in R

[![R-CMD-check](https://github.com/yht/dimvR/actions/workflows/r-check.yml/badge.svg)](https://github.com/yht/dimvR/actions/workflows/r-check.yml)
[![Smoke Benchmark](https://github.com/yht/dimvR/actions/workflows/smoke-benchmark.yml/badge.svg)](https://github.com/yht/dimvR/actions/workflows/smoke-benchmark.yml)

This R package implements a regularized conditional distribution-based approach for missing data imputation,
commonly referred to as **DIMV (Distribution-based Imputation using Conditional Expectation with Regularization)**
[(Nguyen et al., 2023)](https://arxiv.org/abs/2302.00911). The method imputes missing entries by estimating the
conditional expectation of each variable given the others, while incorporating regularization to improve
numerical stability and robustness in high-dimensional or multicollinear settings.

The development of this package is motivated by methodological advances discussed in the paper
*“Explainability of Machine Learning Models under Missing Data”* [(Nguyen et al., 2024)](https://arxiv.org/abs/2407.00411v1),
which highlights the importance of principled imputation methods for maintaining both predictive performance and
model interpretability. The **dimvR** package offers a faithful and reproducible R implementation of the
[DIMV framework](https://github.com/maianhpuco/DIMVImputation), enabling researchers and practitioners to
incorporate regularized conditional imputation into statistical modeling workflows, machine learning pipelines,
and explainability studies.


## Key Features

- **Regularized Conditional Imputation (DIMV):**  
  Implements conditional expectation–based imputation with ridge regularization for enhanced stability in
  high-dimensional or multicollinear data settings.

- **Deterministic and Multiple Imputation:**  
  Supports single deterministic imputation as well as multiple imputation by injecting residual uncertainty to
  better reflect posterior variability.

- **Convergence-Based Iterative Procedure:**  
  Iteratively updates missing values variable-by-variable until convergence, ensuring stable and consistent
  imputations.

- **Feature Selection (Optional):**  
  Adaptive, fixed-threshold, mutual information, and hybrid selection strategies to reduce noise and
  multicollinearity in conditional models.

- **Explainability Pipeline (Optional):**  
  End-to-end evaluation with SHAP-based explainability and HTML report generation (requires optional dependencies).

- **Lightweight Core, Optional Extensions:**  
  Core imputation is implemented in base R; advanced pipelines use optional packages for modeling and visualization.

- **Model-Agnostic Integration:**  
  Output can be used with any downstream statistical or machine learning model, including regression, tree-based
  models, and modern explainability frameworks.

## Engineering Progress (P0)

| Item | Status | Evidence |
|---|---|---|
| CI matrix (Debian/Windows, R release) | In progress | `.github/workflows/r-check.yml` |
| Standard test harness (`testthat`) | Implemented | `tests/testthat.R`, `DESCRIPTION` (`Config/testthat/edition: 3`) |
| Automated smoke benchmark in CI | Implemented | `.github/workflows/smoke-benchmark.yml`, `eval/smoke_benchmark.R` |
| Benchmark artifact publication | Implemented | CI artifact: `ci_smoke_benchmark.csv`, `ci_smoke_summary.md` |
| Progress metrics (test count + coverage) | Implemented | CI artifact: `ci_progress_metrics.csv`, `ci_progress_summary.md` |
| Interim coverage gate | Implemented | `MIN_COVERAGE=40` via `eval/check_coverage_gate.R` in CI workflow |

Latest snapshot (2026-02-09):
- Test files: 8
- Test cases: 25
- Estimated coverage: 59.32%
- Coverage gate status: pass against interim threshold (40%)


## Assumptions & Limitations

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

## Documentation

- Vignette: `vignettes/feature_selection.Rmd`
- Examples: `examples/day1_feature_selection_demo.R`, `examples/benchmark_feature_selection.R`


## To Do

* Method Enhancements
* Explainability & XAI Integration
* Benchmarking & Evaluation Suite
* Production & CRAN Readiness
