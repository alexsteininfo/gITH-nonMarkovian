#!/usr/bin/env bash
# EvoTracer inference stage runner.
#
# Regenerates data/inference/{mutationrate,selection}/<scenario>/{full_trees,
# subsampled_trees}/*.jls by running the 12 driver scripts, then produces the
# 16 figures under figures/3_inference/. Every driver is resumable via
# `isfile(out_path)` checks, so an interrupted run resumes on the next call and
# a forced redo is a file deletion (e.g. `rm -rf data/inference/mutationrate/`).
#
# Usage:
#   bash code/3_inference/run_pipeline.sh              # everything (drivers + figures)
#   bash code/3_inference/run_pipeline.sh drivers      # just data/inference/
#   bash code/3_inference/run_pipeline.sh figures      # just PNGs (needs data)
#   bash code/3_inference/run_pipeline.sh mutrate      # mutation-rate drivers + figures
#   bash code/3_inference/run_pipeline.sh selection    # selection drivers + figures
#
# Each Julia invocation logs to data/inference/logs/<script-stem>.log so a
# multi-hour run in tmux stays readable. `set -o pipefail` propagates Julia
# errors through the `tee`.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

# Julia 1.12.2 is the module default AND the version data/raw/ was serialized
# with (data format v30 — Julia 1.10 tops out at v24 and cannot read it). So
# we stay on 1.12 and mitigate its known parallel-GC corruption bug in
# --gcthreads below.
module load Julia

# Cap compute threads at 8. `-t auto` on this 56-core host spawns 56 threads,
# which both stresses the GC and monopolises the shared HPC. 8 is a good
# middle: driver hot loops scale roughly linearly to ~8 threads and then
# saturate on Serialization + statistic allocation.
THREADS="${JULIA_THREADS:-8}"

# Parallel GC re-enabled after the sister_depth_ratios nperm=999 → 0 fix.
# The original 1.12 GC corruption was triggered by that statistic's
# pathological allocation pattern (~10^11 shuffle-derived StatValue allocs on
# an N16384 shard). With that fixed, allocation pressure is ~1000× lower and
# parallel GC (Julia's default) is safe again. If corruption recurs, revert
# to --gcthreads=1,0.
GC_ARGS=()

STAGE="${1:-all}"

LOG_DIR="$ROOT/data/inference/logs"
mkdir -p "$LOG_DIR"

run_julia() {
    local script="$1"
    local logname
    logname="$(basename "$script" .jl)"
    echo
    echo "==> $(date +%H:%M:%S)  julia --project=. -t $THREADS ${GC_ARGS[*]:-} $script"
    julia --project=. -t "$THREADS" ${GC_ARGS[@]+"${GC_ARGS[@]}"} "$script" 2>&1 | tee "$LOG_DIR/${logname}.log"
}

# Mutation-rate drivers — full first, then subsampled. Within each, neutral
# before sel_{1,2} because the subsampled drivers for sel_{1,2} read matching
# files from data/processed/selection_{1,2}/<timing>/{injection,n_restarts}/,
# not from mutrate outputs — so cross-stage order is not load-bearing here.
MUT_DRIVERS=(
    code/3_inference/01_mutationrate/full_trees/mutation_rate_neutral.jl
    code/3_inference/01_mutationrate/full_trees/mutation_rate_sel1.jl
    code/3_inference/01_mutationrate/full_trees/mutation_rate_sel2.jl
    code/3_inference/01_mutationrate/subsampled_trees/mutation_rate_neutral_subsampled.jl
    code/3_inference/01_mutationrate/subsampled_trees/mutation_rate_sel1_subsampled.jl
    code/3_inference/01_mutationrate/subsampled_trees/mutation_rate_sel2_subsampled.jl
)

SEL_DRIVERS=(
    code/3_inference/02_selection/full_trees/selection_neutral.jl
    code/3_inference/02_selection/full_trees/selection_sel1.jl
    code/3_inference/02_selection/full_trees/selection_sel2.jl
    code/3_inference/02_selection/subsampled_trees/selection_neutral_subsampled.jl
    code/3_inference/02_selection/subsampled_trees/selection_sel1_subsampled.jl
    code/3_inference/02_selection/subsampled_trees/selection_sel2_subsampled.jl
)

# Figures are globbed so this stays in sync as scripts are added or renamed.
MUT_FIGURES=(code/3_inference/01_mutationrate/figures/*.jl)
SEL_FIGURES=(code/3_inference/02_selection/figures/*.jl)

run_all() {
    local -n arr=$1
    for script in "${arr[@]}"; do
        run_julia "$script"
    done
}

case "$STAGE" in
    all)
        run_all MUT_DRIVERS
        run_all SEL_DRIVERS
        run_all MUT_FIGURES
        run_all SEL_FIGURES
        ;;
    drivers)
        run_all MUT_DRIVERS
        run_all SEL_DRIVERS
        ;;
    figures)
        run_all MUT_FIGURES
        run_all SEL_FIGURES
        ;;
    mutrate)
        run_all MUT_DRIVERS
        run_all MUT_FIGURES
        ;;
    selection)
        run_all SEL_DRIVERS
        run_all SEL_FIGURES
        ;;
    *)
        echo "unknown stage: $STAGE" >&2
        echo "usage: $(basename "$0") [all|drivers|figures|mutrate|selection]" >&2
        exit 1
        ;;
esac

echo
echo "run_pipeline.sh: $STAGE complete at $(date +%H:%M:%S)"
