#!/usr/bin/env bash
# plot_phylogeny_examples.sh
#
# Render a small gallery of MEDICC2-vs-truth phylogeny comparisons (paired
# ggtree panels + ape::cophyloplot tanglegram) spanning the same 8 scenarios
# used for the CN-profile input examples: the three neutral timing models,
# death-rate variation, a larger-N run, and one example from each selection
# scenario.
#
# Output: figures/3_inference/03_MEDICC2/phylogenies/{phylogenies_paired,phylogenies_tanglegram}__<slug>.png
#
# Usage:
#   bash code/3_inference/03_MEDICC2/assessment/plot_phylogeny_examples.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
SCRIPT="${REPO_ROOT}/code/3_inference/03_MEDICC2/assessment/plot_phylogenies.R"
OUTDIR="figures/3_inference/03_MEDICC2/phylogenies"  # resolved against REPO_ROOT inside the R script

SIMS=(
    # 3 neutral timing models @ N=1000, d=0.5 (deterministic has no d parameter)
    "neutral/deterministic/neutral_deterministic_N1024_n102/sim1"
    "neutral/gamma/neutral_gamma_N1000_d0.5_k5.0_n100/sim1"
    "neutral/markov/neutral_markov_N1000_d0.5_n100/sim1"
    # Death-rate variation (gamma, N=1000)
    "neutral/gamma/neutral_gamma_N1000_d0.0_k5.0_n100/sim1"
    "neutral/gamma/neutral_gamma_N1000_d0.9_k5.0_n100/sim1"
    # Larger N -> larger MEDICC2 sample size (n=164)
    "neutral/deterministic/neutral_deterministic_N16384_n164/sim1"
    # Selection scenarios
    "selection_1/gamma/sel1_gamma_N1000_d0.5_k5.0_s1.0_n100/sim1"
    "selection_2/gamma/sel2_gamma_N1000_d0.5_k5.0_s0.15_M10.0_n100/sim1"
)

cd "$REPO_ROOT"
for sim in "${SIMS[@]}"; do
    echo "=== $sim ==="
    micromamba run -n R Rscript "$SCRIPT" --sim "$sim" --outdir "$OUTDIR"
done
