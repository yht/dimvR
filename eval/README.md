# CI Smoke Benchmark

This folder contains the script used by GitHub Actions to run a lightweight
benchmark smoke test:

- Script: `eval/smoke_benchmark.R`
- Progress metrics script: `eval/ci_progress_metrics.R`
- Coverage gate script: `eval/check_coverage_gate.R`
- Artifact outputs:
  - `eval/ci_smoke_benchmark.csv`
  - `eval/ci_smoke_summary.md`
  - `eval/ci_progress_metrics.csv`
  - `eval/ci_progress_summary.md`
  - `eval/ci_coverage_gate_summary.md`

The smoke benchmark is intentionally small and fast. It is used to prove that
the benchmarking pipeline is alive and reproducible in CI while larger
multi-dataset benchmarks are still in progress.

Current interim gate: `MIN_COVERAGE=40` (configured in workflow).

Latest checked-in coverage validation (2026-09-14):
- test_files: 14
- test_cases: 57
- coverage_percent: 72.83
- coverage gate: PASS at the interim threshold of 40%

Eight optional/conditional tests were skipped in this environment: the XGBoost/SHAP paths were unavailable, and the missing-backend test was not applicable because `ranger` was installed.

The smoke benchmark artifact is from 2026-09-10; the coverage and progress artifacts were regenerated on 2026-09-14.
