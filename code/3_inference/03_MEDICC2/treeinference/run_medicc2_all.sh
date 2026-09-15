#!/usr/bin/env bash
# Run MEDICC2 on every sim*_cells.tsv file under data/CN_subsampled/.
#
# Input layout (one MEDICC2-ready TSV per simulation):
#   data/CN_subsampled/<scenario>/<model>/<param_dir>/sim<K>_cells.tsv
#
# Output layout (mirrors the input, one directory per simulation):
#   data/MEDICC2/treeinference/<scenario>/<model>/<param_dir>/sim<K>/sim<K>_*
#
# The runner is resumable: a simulation is skipped if
# <output_dir>/<prefix>_final_tree.new already exists. Delete that file (or
# the whole per-sim output directory) to force a rerun.
#
# Run in a dedicated tmux session. N_PARALLEL invocations run simultaneously,
# each using N_CORES cores (~N_PARALLEL * N_CORES total). The HMM burst is the
# only phase that scales with --n-cores; steady-state is ~1 CPU/process, so
# many parallel low-core invocations give better throughput than one high-core
# sequential run.

set -uo pipefail

# Repo root -- climb four dirs from this script's location
# (code/3_inference/03_MEDICC2/treeinference/run_medicc2_all.sh -> repo root).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

INPUT_ROOT="$ROOT/data/CN_subsampled"
OUTPUT_ROOT="$ROOT/data/MEDICC2/treeinference"
LOG_DIR="$OUTPUT_ROOT/logs"

MEDICC2="/srv/home/aste0033/micromamba/envs/python/bin/medicc2"
N_CORES=9
N_PARALLEL=4

mkdir -p "$OUTPUT_ROOT" "$LOG_DIR"

if [[ ! -x "$MEDICC2" ]]; then
    echo "ERROR: medicc2 not found at $MEDICC2" >&2
    exit 1
fi
if [[ ! -d "$INPUT_ROOT" ]]; then
    echo "ERROR: input root not found: $INPUT_ROOT" >&2
    exit 1
fi

# Collect all *_cells.tsv inputs in a stable order.
mapfile -t INPUTS < <(find "$INPUT_ROOT" -type f -name '*_cells.tsv' | sort)
TOTAL=${#INPUTS[@]}
if [[ $TOTAL -eq 0 ]]; then
    echo "ERROR: no *_cells.tsv files found under $INPUT_ROOT" >&2
    exit 1
fi

echo "MEDICC2 full sweep"
echo "  inputs:       $TOTAL"
echo "  parallel:     $N_PARALLEL slots"
echo "  cores/slot:   $N_CORES"
echo "  input root:   $INPUT_ROOT"
echo "  output root:  $OUTPUT_ROOT"
echo

run_one() {
    local input="$1"
    local fname prefix rel_param_dir out_dir sentinel log_file t0 t1 dt

    fname="$(basename "$input")"
    prefix="${fname%_cells.tsv}"
    rel_param_dir="${input#$INPUT_ROOT/}"
    rel_param_dir="$(dirname "$rel_param_dir")"
    out_dir="$OUTPUT_ROOT/$rel_param_dir/$prefix"
    sentinel="$out_dir/${prefix}_final_tree.new"

    if [[ -f "$sentinel" ]]; then
        printf '[skip  %s] %s\n' "$(date +%H:%M:%S)" "$rel_param_dir/$prefix"
        return 0
    fi

    mkdir -p "$out_dir"
    log_file="$LOG_DIR/${rel_param_dir//\//__}__${prefix}.log"
    mkdir -p "$(dirname "$log_file")"

    printf '[start %s] %s\n' "$(date +%H:%M:%S)" "$rel_param_dir/$prefix"
    t0=$(date +%s)
    if "$MEDICC2" \
            "$input" \
            "$out_dir" \
            --input-type t \
            --prefix "$prefix" \
            --n-cores "$N_CORES" \
            --events \
            --plot none \
            >"$log_file" 2>&1; then
        t1=$(date +%s); dt=$((t1 - t0))
        if [[ -f "$sentinel" ]]; then
            printf '[ok    %s] %s  (%ds)\n' "$(date +%H:%M:%S)" "$rel_param_dir/$prefix" "$dt"
        else
            printf '[fail  %s] %s  (%ds, no sentinel) see %s\n' \
                "$(date +%H:%M:%S)" "$rel_param_dir/$prefix" "$dt" "$log_file"
        fi
    else
        t1=$(date +%s); dt=$((t1 - t0))
        printf '[fail  %s] %s  (%ds, exit non-zero) see %s\n' \
            "$(date +%H:%M:%S)" "$rel_param_dir/$prefix" "$dt" "$log_file"
    fi
    return 0
}
export -f run_one
export MEDICC2 N_CORES INPUT_ROOT OUTPUT_ROOT LOG_DIR

RUN_START=$(date +%s)
printf '%s\0' "${INPUTS[@]}" \
    | xargs -0 -n1 -P "$N_PARALLEL" -I {} bash -c 'run_one "$@"' _ {}
RUN_END=$(date +%s)

echo
echo "Elapsed: $((RUN_END - RUN_START))s"
