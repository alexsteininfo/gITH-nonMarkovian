#!/usr/bin/env bash
# Reduced-scope EvoTracer inference-stage runner.
#
# Runs every driver + every figure in run_pipeline.sh, but processes at most
# $INFERENCE_MAX_SIMS simulations per shard. Sims past the cap remain `nothing`
# in the output vector (figure code already skips these).
#
# Output paths are IDENTICAL to run_pipeline.sh. So after a reduced run the
# `.jls` files under data/inference/ look "complete" to the resumability check;
# the full pipeline would skip them silently. **Clean data/inference/ before
# switching from reduced to full**:
#
#     rm -rf data/inference
#     bash code/3_inference/run_pipeline.sh
#
# Usage:
#   bash code/3_inference/run_pipeline_reduced.sh                  # sims=20, all stages
#   INFERENCE_MAX_SIMS=5 bash .../run_pipeline_reduced.sh          # even smaller
#   bash .../run_pipeline_reduced.sh drivers                       # just the .jls
#   bash .../run_pipeline_reduced.sh figures                       # just the PNGs
#   bash .../run_pipeline_reduced.sh mutrate                       # mut-rate half
#   bash .../run_pipeline_reduced.sh selection                     # selection half
#
# Default cap is 20 sims — matches the ballpark you'd get from a rough eyeball
# of a small parameter sweep, small enough to finish across the full grid in
# under an hour, big enough that the figures are not degenerate.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

export INFERENCE_MAX_SIMS="${INFERENCE_MAX_SIMS:-20}"

echo "run_pipeline_reduced.sh: INFERENCE_MAX_SIMS=$INFERENCE_MAX_SIMS"
echo "  (all statistics run; only the first $INFERENCE_MAX_SIMS sims per shard)"
echo "  output paths are the same as run_pipeline.sh — clean data/inference/ before running full."
echo

exec bash "$ROOT/code/3_inference/run_pipeline.sh" "$@"
