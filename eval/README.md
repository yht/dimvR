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

Latest local baseline (2026-02-09):
- test_files: 8
- test_cases: 25
- coverage_percent: 59.32
