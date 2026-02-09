# CI Smoke Benchmark

This folder contains the script used by GitHub Actions to run a lightweight
benchmark smoke test:

- Script: `eval/smoke_benchmark.R`
- Artifact outputs:
  - `eval/ci_smoke_benchmark.csv`
  - `eval/ci_smoke_summary.md`

The smoke benchmark is intentionally small and fast. It is used to prove that
the benchmarking pipeline is alive and reproducible in CI while larger
multi-dataset benchmarks are still in progress.
