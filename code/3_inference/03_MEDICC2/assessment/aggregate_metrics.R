# aggregate_metrics.R
#
# For every MEDICC2 sim on disk, compute a small set of summary statistics
# and cache to CSV. Consumed by plot_aggregate_metrics.R.
#
# Statistics per sim:
#   - truth_bl_mean : mean events per collapsed truth edge
#   - med_bl_mean   : mean events per MEDICC2 edge
#   - truth_EM      : mean per-cell mutational burden (root->leaf events, truth)
#   - med_EM        : mean per-cell mutational burden (MEDICC2, diploid dropped)
#   - coph_slope    : slope of MEDICC2 cophenetic pairwise distances regressed
#                     on truth cophenetic pairwise distances
#   - rf_raw        : Robinson-Foulds distance (unnormalised)
#   - rf_norm       : Robinson-Foulds distance (normalised to [0, 1])
#
# Definitions match plot_saturation.R and plot_topology_metrics.R exactly, so
# per-sim aggregates equal what the per-sim panels show.
#
# Cache: data/MEDICC/benchmark/aggregate_metrics.csv (created if missing).
# Idempotent: sims already present in the CSV are skipped unless --force.
# Backfill: if the CSV exists but lacks any of the metric columns above (e.g.
# a new column was added), rows missing that value are recomputed in place
# rather than requiring --force.
#
# Usage:
#   Rscript code/3_inference/03_MEDICC2/assessment/aggregate_metrics.R [--force]

suppressPackageStartupMessages({
    library(dplyr); library(readr); library(tibble); library(ape); library(phangorn)
})

# Metric columns produced by compute_metrics_for_sim(). Kept as a single list
# so the backfill logic and the failure-path NA row can stay in sync.
METRIC_COLS <- c("truth_bl_mean", "med_bl_mean", "truth_EM", "med_EM",
                 "coph_slope", "rf_raw", "rf_norm")

args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE))
if (length(this_file) != 1L) stop("cannot resolve --file")
ROOT <- normalizePath(file.path(dirname(this_file), "..", "..", "..", ".."))

source(file.path(ROOT, "code", "X_helpers", "medicc2_io.R"))

user_args <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% user_args

CSV_PATH <- file.path(ROOT, "data", "MEDICC", "benchmark", "aggregate_metrics.csv")
dir.create(dirname(CSV_PATH), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Sim discovery
# ---------------------------------------------------------------------------
MED_ROOT <- file.path(ROOT, "data", "MEDICC2", "treeinference")
tree_files <- list.files(MED_ROOT, pattern = "^sim.*_final_tree\\.new$",
                         recursive = TRUE, full.names = TRUE)
sim_rels <- unique(dirname(sub(paste0("^", MED_ROOT, "/"), "", tree_files)))
sim_rels <- sort(sim_rels)
message(sprintf("Discovered %d sims under %s", length(sim_rels), MED_ROOT))

# ---------------------------------------------------------------------------
# Parameter parser: pulls N/d/k/s/M/n out of the parameter-dir name
# ---------------------------------------------------------------------------
parse_sim_params <- function(sim_rel) {
    parts <- strsplit(sim_rel, "/", fixed = TRUE)[[1]]
    stopifnot(length(parts) >= 4)
    param_dir <- parts[3]
    ex_num <- function(pat) {
        m <- regmatches(param_dir, regexec(pat, param_dir))[[1]]
        if (length(m) < 2) NA_real_ else as.numeric(m[2])
    }
    tibble(
        sim_rel  = sim_rel,
        scenario = parts[1],
        timing   = parts[2],
        sim_id   = parts[4],
        N     = as.integer(ex_num("_N([0-9]+)")),
        d     = ex_num("_d([0-9.]+)"),
        k     = ex_num("_k([0-9.]+)"),
        s     = ex_num("_s([0-9.]+)"),
        M_inj = ex_num("_M([0-9.]+)"),
        n     = as.integer(ex_num("_n([0-9]+)"))
    )
}

# ---------------------------------------------------------------------------
# Per-sim metric computation (definitions mirror plot_saturation.R and
# plot_topology_metrics.R, so aggregates equal what the per-sim panels show)
# ---------------------------------------------------------------------------
compute_metrics_for_sim <- function(sim_rel) {
    truth_dir <- file.path(ROOT, "data", "CN_subsampled", dirname(sim_rel))
    med_dir   <- file.path(ROOT, "data", "MEDICC2", "treeinference", sim_rel)

    truth_tree_raw    <- read_truth_tree(truth_dir)
    truth_events      <- read_truth_events(truth_dir)
    medicc2_tree_full <- read_medicc2_tree(med_dir)
    bl_df             <- read_medicc2_branch_lengths(med_dir)

    collapsed  <- collapse_unary(truth_tree_raw)
    truth_tree <- collapsed$tree
    ec_truth   <- edge_event_counts_truth(collapsed, truth_events)
    # Assign truth branch lengths in "events" so cophenetic uses the same unit
    # as MEDICC2's own branch lengths (matches plot_topology_metrics.R).
    truth_tree$edge.length <- ec_truth$events

    bt <- per_cell_burden_truth(collapsed, truth_events)

    medicc2_tree <- if ("diploid" %in% medicc2_tree_full$tip.label)
        ape::drop.tip(medicc2_tree_full, "diploid") else medicc2_tree_full
    bm <- per_cell_burden_medicc2(medicc2_tree_full, bl_df) %>%
        filter(cell != "diploid")

    truth_bl_mean <- mean(as.numeric(ec_truth$events))
    med_bl_mean   <- mean(as.numeric(bl_df$length))
    truth_EM      <- mean(as.numeric(bt$M))
    med_EM        <- mean(as.numeric(bm$M))

    common <- intersect(truth_tree$tip.label, medicc2_tree$tip.label)
    if (length(common) < 0.9 * length(truth_tree$tip.label)) {
        coph_slope <- NA_real_
        rf_raw <- NA_real_
        rf_norm <- NA_real_
    } else {
        truth_c <- if (length(common) < length(truth_tree$tip.label))
            ape::drop.tip(truth_tree, setdiff(truth_tree$tip.label, common)) else truth_tree
        med_c   <- if (length(common) < length(medicc2_tree$tip.label))
            ape::drop.tip(medicc2_tree, setdiff(medicc2_tree$tip.label, common)) else medicc2_tree
        co_t <- ape::cophenetic.phylo(truth_c)[common, common]
        co_m <- ape::cophenetic.phylo(med_c)[common, common]
        uti  <- upper.tri(co_t)
        fit  <- lm(co_m[uti] ~ co_t[uti])
        coph_slope <- unname(coef(fit)[2])

        # Robinson-Foulds (topology-only). Matches plot_topology_metrics.R.
        rf_raw  <- tryCatch(as.numeric(phangorn::RF.dist(truth_c, med_c, normalize = FALSE)),
                            error = function(e) NA_real_)
        rf_norm <- tryCatch(as.numeric(phangorn::RF.dist(truth_c, med_c, normalize = TRUE)),
                            error = function(e) NA_real_)
    }

    tibble(
        truth_bl_mean = truth_bl_mean,
        med_bl_mean   = med_bl_mean,
        truth_EM      = truth_EM,
        med_EM        = med_EM,
        coph_slope    = coph_slope,
        rf_raw        = rf_raw,
        rf_norm       = rf_norm
    )
}

# ---------------------------------------------------------------------------
# Load existing CSV, decide what still needs computing
#
# A sim needs computing if it is not in the CSV at all, OR if any of the
# METRIC_COLS is missing / NA for that row (backfill for newly added metrics).
# ---------------------------------------------------------------------------
existing <- if (file.exists(CSV_PATH) && !force) {
    suppressMessages(read_csv(CSV_PATH, show_col_types = FALSE))
} else {
    tibble()
}

needs_backfill <- if (nrow(existing) > 0L) {
    missing_cols <- setdiff(METRIC_COLS, names(existing))
    if (length(missing_cols) > 0L) {
        message(sprintf("CSV lacks columns: %s -- backfilling all %d cached rows",
                        paste(missing_cols, collapse = ", "), nrow(existing)))
        existing$sim_rel
    } else {
        na_rows <- rowSums(is.na(existing[, METRIC_COLS])) > 0L
        existing$sim_rel[na_rows]
    }
} else {
    character()
}

already <- if (nrow(existing) > 0L) setdiff(existing$sim_rel, needs_backfill) else character()
todo <- sort(unique(c(setdiff(sim_rels, already), intersect(sim_rels, needs_backfill))))
message(sprintf("Already cached: %d, to (re)compute: %d (of which %d backfill), force=%s",
                length(already), length(todo),
                length(intersect(todo, needs_backfill)), force))

# ---------------------------------------------------------------------------
# Compute
# ---------------------------------------------------------------------------
new_rows <- vector("list", length(todo))
t0 <- Sys.time()
for (i in seq_along(todo)) {
    sim_rel <- todo[i]
    params <- parse_sim_params(sim_rel)
    metrics <- tryCatch(compute_metrics_for_sim(sim_rel),
                        error = function(e) {
                            message(sprintf("  [%3d/%d] FAILED %s: %s",
                                            i, length(todo), sim_rel, e$message))
                            as_tibble(setNames(
                                replicate(length(METRIC_COLS), NA_real_, simplify = FALSE),
                                METRIC_COLS
                            ))
                        })
    new_rows[[i]] <- bind_cols(params, metrics)
    if (i %% 10L == 0L || i == length(todo)) {
        dt <- as.numeric(Sys.time() - t0, units = "secs")
        message(sprintf("  [%3d/%d] %.1fs elapsed (%.2fs/sim avg)",
                        i, length(todo), dt, dt / i))
    }
}

# ---------------------------------------------------------------------------
# Merge + atomic write
#
# Backfilled sims appear both in `existing` and in `new_df`; drop the stale
# rows from `existing` so the recomputed values win rather than duplicate.
# ---------------------------------------------------------------------------
new_df <- if (length(new_rows)) bind_rows(new_rows) else tibble()
combined <- if (force || nrow(existing) == 0L) {
    new_df
} else {
    kept <- existing %>% filter(!(sim_rel %in% new_df$sim_rel))
    bind_rows(kept, new_df)
}
combined <- combined %>% arrange(scenario, timing, N, d, s, M_inj)

tmp_path <- paste0(CSV_PATH, ".tmp")
write_csv(combined, tmp_path)
file.rename(tmp_path, CSV_PATH)
message(sprintf("wrote: %s (%d rows)", CSV_PATH, nrow(combined)))
