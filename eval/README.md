# CI Smoke Benchmark

This folder is the active checked-in evaluation and progress baseline for the
package source tree. It contains the scripts used by GitHub Actions to run a
lightweight benchmark smoke test, summarize test/progress metrics, and enforce
the interim coverage gate:

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

Authoritative checked-in snapshot (2026-04-06):
- test_files: 11
- test_cases: 34
- coverage_percent: 76.87

Status notes:
- `dimvExplainR/eval/` is the active reporting baseline for package work.
- The root-level `../eval/` folder is archival/stale unless it is explicitly regenerated.
- Refresh the artifacts in this folder after meaningful test or coverage changes, then sync the resulting snapshot into `PROGRESS.md` and any active package-facing docs.
