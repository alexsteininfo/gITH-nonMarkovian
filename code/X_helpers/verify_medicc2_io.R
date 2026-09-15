# verify_medicc2_io.R
#
# Smoke test for medicc2_io.R helpers against the target sim.
# Usage:
#   micromamba run -n R Rscript code/X_helpers/verify_medicc2_io.R

suppressPackageStartupMessages({
    library(dplyr)
})

# Resolve repo root from --file= (matches plot_trees.R idiom)
args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
if (length(this_file) != 1L) stop("verify_medicc2_io.R: cannot resolve --file")
ROOT <- normalizePath(file.path(dirname(this_file), "..", ".."))

source(file.path(ROOT, "code", "X_helpers", "medicc2_io.R"))

SIM <- "neutral/gamma/neutral_gamma_N1000_d0.5_k5.0_n100/sim1"
# Truth files live directly in the param dir (no sim1/ subdirectory):
#   data/CN_subsampled/neutral/gamma/neutral_gamma_N1000_d0.5_k5.0_n100/sim1_*.{tsv,nwk}
TRUTH_DIR <- file.path(ROOT, "data", "CN_subsampled", dirname(SIM))
# MEDICC2 output has the sim1/ level present:
#   data/MEDICC2/treeinference/neutral/gamma/.../sim1/
MEDICC2_DIR <- file.path(ROOT, "data", "MEDICC2",
                          "treeinference", dirname(SIM), basename(SIM))

message("truth dir:   ", TRUTH_DIR)
message("medicc2 dir: ", MEDICC2_DIR)
stopifnot(dir.exists(TRUTH_DIR), dir.exists(MEDICC2_DIR))

# --- Loaders ---
truth_tree_raw <- read_truth_tree(TRUTH_DIR)
medicc2_tree   <- read_medicc2_tree(MEDICC2_DIR)
truth_events   <- read_truth_events(TRUTH_DIR)
bl_df          <- read_medicc2_branch_lengths(MEDICC2_DIR)

stopifnot(inherits(truth_tree_raw, "phylo"),
          length(truth_tree_raw$tip.label) > 0L)
stopifnot(inherits(medicc2_tree, "phylo"),
          length(medicc2_tree$tip.label) > 0L)
stopifnot(nrow(truth_events) > 0L, nrow(bl_df) > 0L)

# --- collapse_unary sanity ---
collapsed <- collapse_unary(truth_tree_raw)
stopifnot(inherits(collapsed$tree, "phylo"))
stopifnot(length(collapsed$chain_nodes) == nrow(collapsed$tree$edge))
stopifnot(all(collapsed$divisions_collapsed >= 1L))
message(sprintf("collapsed: %d tips (was %d), %d edges",
                length(collapsed$tree$tip.label),
                length(truth_tree_raw$tip.label),
                nrow(collapsed$tree$edge)))

# The collapsed truth tree's tips should be a subset of the MEDICC2 tree's
# tips (MEDICC2 adds `diploid`). Both should share the 100 sampled cells.
truth_tips <- collapsed$tree$tip.label
medicc2_tips <- medicc2_tree$tip.label
common <- intersect(truth_tips, medicc2_tips)
message(sprintf("common tips: %d (truth=%d, medicc2=%d)",
                length(common), length(truth_tips), length(medicc2_tips)))
if (length(common) < 0.9 * length(truth_tips)) {
    stop("verify: <90% of truth tips are present in the MEDICC2 tree")
}

# --- Burden functions ---
bt <- per_cell_burden_truth(collapsed, truth_events)
bm <- per_cell_burden_medicc2(medicc2_tree, bl_df)
stopifnot(nrow(bt) == length(truth_tips), nrow(bm) == length(medicc2_tips))
stopifnot(all(bt$M >= 0L), all(bm$M >= 0))
message(sprintf("truth   M: min=%d  median=%.1f  max=%d  E=%.2f  Var=%.2f",
                min(bt$M), median(bt$M), max(bt$M), mean(bt$M), var(bt$M)))
message(sprintf("medicc2 M: min=%.1f  median=%.1f  max=%.1f  E=%.2f  Var=%.2f",
                min(bm$M), median(bm$M), max(bm$M), mean(bm$M), var(bm$M)))

# --- Edge event counts ---
ec <- edge_event_counts_truth(collapsed, truth_events)
stopifnot(nrow(ec) == nrow(collapsed$tree$edge))
message(sprintf("edge events (truth): total=%d  max=%d",
                sum(ec$events), max(ec$events)))
message(sprintf("edge lengths (medicc2): total=%.0f  max=%.0f",
                sum(bl_df$length), max(bl_df$length)))

# --- Stratified sample ---
picks <- stratified_cell_sample(bt, n_stratified = 4L, n_random = 2L, seed = 42L)
stopifnot(length(picks) == 6L, length(unique(picks)) == 6L, all(picks %in% bt$cell))
message("stratified_cell_sample picks (seeded): ", paste(picks, collapse = ", "))

message("\nOK: medicc2_io helpers pass smoke tests")
