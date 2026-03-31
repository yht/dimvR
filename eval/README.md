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

Latest checked-in CI baseline (2026-03-30):
- test_files: 11
- test_cases: 30
- coverage_percent: NA

Working-tree review during documentation pass (2026-03-31):
- test_files: 10
- test_cases: 29
- coverage_percent: not re-measured

These numbers indicate that the checked-in progress artifacts should be regenerated so the repository metrics match the current tree.
