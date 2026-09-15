#!/usr/bin/env bash
# Curated, 3-way parallel MEDICC2 runner.
#
# Runs a hand-picked subset of param dirs under data/CN_subsampled/, three
# MEDICC2 invocations at a time with --n-cores 16 each (~48 of 56 cores).
# Exists so the sel1 subset needed for the 16:00 supervisor slot can complete
# in an afternoon rather than the ~days a sequential full-sweep would take.
#
# Layout, sentinel, and log conventions are the same as run_medicc2_all.sh, so
# outputs from this runner are indistinguishable from the sweep's outputs, and
# the sentinel skip means already-done items (all completed neutral runs) are
# re-checked and skipped instantly.
#
# Run in a dedicated tmux session; do not launch multiple copies in parallel.

set -uo pipefail

# Repo root -- climb four dirs from this script's location
# (code/3_inference/03_MEDICC2/treeinference/run_medicc2_curated.sh -> repo root).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

INPUT_ROOT="$ROOT/data/CN_subsampled"
OUTPUT_ROOT="$ROOT/data/MEDICC2/treeinference"
LOG_DIR="$OUTPUT_ROOT/logs"

MEDICC2="/srv/home/aste0033/micromamba/envs/python/bin/medicc2"
N_CORES=16
N_PARALLEL=3

# Whitelist of param dirs relative to INPUT_ROOT.
# Only sim1 currently exists per param dir; the runner picks up sim*_cells.tsv.
PARAM_DIRS=(
    selection_1/gamma/sel1_gamma_N1000_d0.0_k5.0_s0.5_n100
    selection_1/gamma/sel1_gamma_N1000_d0.0_k5.0_s1.0_n100
    selection_1/gamma/sel1_gamma_N1000_d0.0_k5.0_s1.5_n100
    selection_1/gamma/sel1_gamma_N1000_d0.0_k5.0_s2.0_n100
    selection_1/gamma/sel1_gamma_N1000_d0.5_k5.0_s1.0_n100
    selection_1/gamma/sel1_gamma_N1000_d0.5_k5.0_s2.0_n100
)

mkdir -p "$OUTPUT_ROOT" "$LOG_DIR"

if [[ ! -x "$MEDICC2" ]]; then
    echo "ERROR: medicc2 not found at $MEDICC2" >&2
    exit 1
fi
if [[ ! -d "$INPUT_ROOT" ]]; then
    echo "ERROR: input root not found: $INPUT_ROOT" >&2
    exit 1
fi

# Resolve whitelist -> actual input files, failing loudly if a param dir has
# gone missing or has no cells.tsv.
INPUTS=()
for rel in "${PARAM_DIRS[@]}"; do
    dir="$INPUT_ROOT/$rel"
    if [[ ! -d "$dir" ]]; then
        echo "ERROR: whitelist entry missing on disk: $dir" >&2
        exit 1
    fi
    matches=()
    while IFS= read -r -d '' f; do
        matches+=("$f")
    done < <(find "$dir" -maxdepth 1 -type f -name '*_cells.tsv' -print0 | sort -z)
    if [[ ${#matches[@]} -eq 0 ]]; then
        echo "ERROR: no *_cells.tsv under $dir" >&2
        exit 1
    fi
    INPUTS+=("${matches[@]}")
done

TOTAL=${#INPUTS[@]}
echo "Curated MEDICC2 sweep"
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
        printf '[skip %s] %s\n' "$(date +%H:%M:%S)" "$rel_param_dir/$prefix"
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
