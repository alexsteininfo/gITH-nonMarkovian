#!/usr/bin/env bash
# plot_parsimony_saturation_examples.sh
#
# Render parsimony-saturation diagnostics + topology-metric figures for the
# same 8 scenarios used by the CN-profile and phylogeny galleries.
#
# Output: figures/3_inference/03_MEDICC2/parsimony_saturation/{saturation,topology_metrics}__<slug>.png
#
# Usage:
#   bash code/3_inference/03_MEDICC2/assessment/plot_parsimony_saturation_examples.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
SCRIPT_DIR="${REPO_ROOT}/code/3_inference/03_MEDICC2/assessment"
OUTDIR="figures/3_inference/03_MEDICC2/parsimony_saturation"  # resolved against REPO_ROOT inside the R scripts

SIMS=(
    "neutral/deterministic/neutral_deterministic_N1024_n102/sim1"
    "neutral/gamma/neutral_gamma_N1000_d0.5_k5.0_n100/sim1"
    "neutral/markov/neutral_markov_N1000_d0.5_n100/sim1"
    "neutral/gamma/neutral_gamma_N1000_d0.0_k5.0_n100/sim1"
    "neutral/gamma/neutral_gamma_N1000_d0.9_k5.0_n100/sim1"
    "neutral/deterministic/neutral_deterministic_N16384_n164/sim1"
    "selection_1/gamma/sel1_gamma_N1000_d0.5_k5.0_s1.0_n100/sim1"
    "selection_2/gamma/sel2_gamma_N1000_d0.5_k5.0_s0.15_M10.0_n100/sim1"
)

cd "$REPO_ROOT"
for sim in "${SIMS[@]}"; do
    echo "=== $sim ==="
    micromamba run -n R Rscript "${SCRIPT_DIR}/plot_saturation.R"       --sim "$sim" --outdir "$OUTDIR"
    micromamba run -n R Rscript "${SCRIPT_DIR}/plot_topology_metrics.R" --sim "$sim" --outdir "$OUTDIR"
done
