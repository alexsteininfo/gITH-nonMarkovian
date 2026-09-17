#!/usr/bin/env bash
# run_aggregate_metrics.sh
#
# Driver for MEDICC2 assessment aggregate figures:
#
#   Parsimony-saturation + phylogeny topology (single scalars per sim):
#     1. aggregate_metrics.R      -- compute + cache per-sim stats to CSV
#        (idempotent; --force to recompute everything; backfills new metric
#        columns automatically).
#     2. plot_aggregate_metrics.R -- read CSV, emit figures under
#        figures/3_inference/03_MEDICC2/{parsimony_saturation,phylogenies}/aggregate/
#
#   Sample-MRCA CN reconstruction:
#     3. mrca_reconstruction_metrics.R -- per-sim bin-level classification of
#        truth vs MEDICC2 at the sample MRCA. Cached to
#        data/MEDICC2/benchmark/mrca_reconstruction.csv. Idempotent;
#        --force to recompute; --cores N for parallelism across sims.
#     4. plot_mrca_reconstruction.R    -- figures under
#        figures/3_inference/03_MEDICC2/MRCA_test/:
#          mrca_classification.png     (truth vs MEDICC2 confusion matrix)
#          mrca_hallucinations.png     (severity of MEDICC2 hallucinations when
#                                       the truth MRCA is diploid)
#          mrca_altered_scatter.png    (bin-level agreement when the truth
#                                       MRCA is non-diploid)
#
# The older `cn_ancestor_metrics.R` + `plot_cn_ancestor.R` pair is
# superseded: it hard-coded `dip = 2L` while both truth and MEDICC2 store CN
# per-haplotype (diploid = 1), so its baseline-relative panels were comparing
# against a nonsensical all-CN=2 reference.
#
# Arguments are forwarded to the two *_metrics.R scripts. Plot scripts take
# no args.
#
# Usage:
#   bash code/3_inference/03_MEDICC2/assessment/run_aggregate_metrics.sh [--force] [--cores N]

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)"
SCRIPT_DIR="${REPO_ROOT}/code/3_inference/03_MEDICC2/assessment"

cd "$REPO_ROOT"
micromamba run -n R Rscript "${SCRIPT_DIR}/aggregate_metrics.R"           "$@"
micromamba run -n R Rscript "${SCRIPT_DIR}/plot_aggregate_metrics.R"
micromamba run -n R Rscript "${SCRIPT_DIR}/mrca_reconstruction_metrics.R" "$@"
micromamba run -n R Rscript "${SCRIPT_DIR}/plot_mrca_reconstruction.R"
