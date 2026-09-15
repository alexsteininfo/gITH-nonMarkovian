#!/usr/bin/env bash
# run_aggregate_metrics.sh
#
# Two-step driver for the aggregate parsimony-saturation summary figures:
#   1. aggregate_metrics.R      -- compute + cache per-sim stats to CSV
#      (idempotent: skips sims already in data/MEDICC/benchmark/aggregate_metrics.csv;
#       pass --force to recompute everything)
#   2. plot_aggregate_metrics.R -- read CSV, emit 3 figures under
#      figures/3_inference/03_MEDICC2/parsimony_saturation/aggregate/
#
# Usage:
#   bash code/3_inference/03_MEDICC2/assessment/run_aggregate_metrics.sh [--force]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
SCRIPT_DIR="${REPO_ROOT}/code/3_inference/03_MEDICC2/assessment"

cd "$REPO_ROOT"
micromamba run -n R Rscript "${SCRIPT_DIR}/aggregate_metrics.R"      "$@"
micromamba run -n R Rscript "${SCRIPT_DIR}/plot_aggregate_metrics.R"
