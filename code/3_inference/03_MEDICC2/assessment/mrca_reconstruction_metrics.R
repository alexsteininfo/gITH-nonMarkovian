# mrca_reconstruction_metrics.R
#
# For every MEDICC2 sim on disk, compare the CN profile MEDICC2 reconstructs at
# the sample-MRCA (root of the tree after dropping the diploid outgroup) against
# the truth CN profile at the MRCA of the sampled cells in the (collapsed) truth
# tree.
#
# The truth format stores per-haplotype CN (one row per (chrom, haplotype)
# segment) and MEDICC2 stores per-haplotype cn_a / cn_b. Both use the convention
# "diploid = 1 per haplotype", so we set dip = 1L throughout. (An earlier
# analysis used dip = 2L, which silently compared everything against an
# all-CN=2 baseline and made every sim look near-perfect regardless of what
# MEDICC2 actually reconstructed.)
#
# For each sim we build two integer vectors of length 2 * n_bins on MEDICC2's
# 1 Mb bin grid: cn_a bins followed by cn_b bins. Every position in the vector
# is one "haplotype-bin". A haplotype-bin is "altered" iff its value differs
# from 1.
#
# Per-sim outputs written to data/MEDICC2/benchmark/mrca_reconstruction.csv,
# one row per sim, with columns:
#   sim_rel, scenario, timing, sim_id, N, d, k, s, M_inj, n
#   truth_mrca_name, medicc_mrca_name, n_hbins   (= 2*n_bins)
#   truth-side counts and CN-unit totals:
#     truth_altered_hbins, truth_altered_units, truth_is_diploid (logical)
#   MEDICC2-side counts and CN-unit totals:
#     med_altered_hbins, med_altered_units, med_is_diploid (logical)
#   overall disagreement:
#     disagree_hbins, disagree_units, exact_match (logical: all bins identical)
#   bin-level confusion (each hbin classified as truth-altered vs MEDICC2-
#   altered):
#     tp_hbins  -- both altered      (correctly flagged as altered)
#     fp_hbins  -- only MEDICC2      (hallucination)
#     fn_hbins  -- only truth        (missed alteration)
#     tn_hbins  -- both diploid      (correctly flagged as diploid)
#   hallucination decomposition (truth diploid at these hbins):
#     halluc_units       -- sum |med - 1| over truth-diploid hbins
#     halluc_gain_units  -- sum (med - 1) over hbins where med > 1
#     halluc_loss_units  -- sum (1 - med) over hbins where med < 1
#   missed decomposition (MEDICC2 diploid at these hbins):
#     missed_units       -- sum |truth - 1| over med-diploid hbins where truth != 1
#   miscall decomposition (both altered, but differ):
#     miscall_hbins, miscall_units
#
# Usage:
#   Rscript code/3_inference/03_MEDICC2/assessment/mrca_reconstruction_metrics.R \
#       [--force] [--cores N] [--sim <sim_rel>]

suppressPackageStartupMessages({
    library(dplyr); library(readr); library(tibble)
    library(ape); library(data.table); library(parallel)
})

args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE))
if (length(this_file) != 1L) stop("cannot resolve --file")
ROOT <- normalizePath(file.path(dirname(this_file), "..", "..", "..", ".."))

source(file.path(ROOT, "code", "X_helpers", "medicc2_io.R"))

# --- CLI --------------------------------------------------------------------
user_args <- commandArgs(trailingOnly = TRUE)
force <- FALSE; cores <- 1L; sim_filter <- NULL
i <- 1L
while (i <= length(user_args)) {
    a <- user_args[i]
    if      (a == "--force") { force <- TRUE; i <- i + 1L }
    else if (a == "--cores") { cores <- as.integer(user_args[i + 1L]); i <- i + 2L }
    else if (a == "--sim")   { sim_filter <- user_args[i + 1L]; i <- i + 2L }
    else stop("unknown arg: ", a)
}
message(sprintf("force=%s cores=%d sim_filter=%s",
                force, cores, if (is.null(sim_filter)) "<all>" else sim_filter))

OUT_CSV <- file.path(ROOT, "data", "MEDICC2", "benchmark",
                     "mrca_reconstruction.csv")
dir.create(dirname(OUT_CSV), showWarnings = FALSE, recursive = TRUE)

MED_ROOT <- file.path(ROOT, "data", "MEDICC2", "treeinference")

# --- Sim discovery ----------------------------------------------------------
tree_files <- list.files(MED_ROOT, pattern = "^sim.*_final_tree\\.new$",
                         recursive = TRUE, full.names = TRUE)
sim_rels <- sort(unique(dirname(sub(paste0("^", MED_ROOT, "/"), "", tree_files))))
if (!is.null(sim_filter)) sim_rels <- intersect(sim_rels, sim_filter)
message(sprintf("Discovered %d sims", length(sim_rels)))

parse_sim_params <- function(sim_rel) {
    parts <- strsplit(sim_rel, "/", fixed = TRUE)[[1]]
    stopifnot(length(parts) >= 4L)
    param_dir <- parts[3]
    ex_num <- function(pat) {
        m <- regmatches(param_dir, regexec(pat, param_dir))[[1]]
        if (length(m) < 2L) NA_real_ else as.numeric(m[2])
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

# --- CN vector builders (MEDICC2 grid) --------------------------------------
build_bin_grid <- function(med_dt) {
    grid <- unique(med_dt[, .(chrom, start, end)])
    setorder(grid, chrom, start)
    grid[, bin_id := seq_len(.N)]
    grid[]
}

# MEDICC2 CN vector for a single node (cn_a then cn_b).
medicc_cn_vec <- function(med_dt, grid, node_name) {
    n_bins <- nrow(grid)
    v <- integer(2L * n_bins) + 1L      # default to diploid (1 per hap)
    sub <- med_dt[sample_id == node_name]
    if (nrow(sub) == 0L) return(NULL)
    sub <- merge(sub, grid, by = c("chrom", "start", "end"), sort = FALSE)
    v[sub$bin_id]           <- sub$cn_a
    v[sub$bin_id + n_bins]  <- sub$cn_b
    v
}

# Truth CN vector for a single node, projected onto MEDICC2 grid via interval
# overlap. When a MEDICC2 bin overlaps several truth segments with different
# CN (which happens whenever a truth-segment boundary falls inside a bin),
# we assign the CN of the truth segment covering the LARGEST fraction of the
# bin -- the dominant state at this resolution. Doing anything else (e.g.
# last-assigned-wins) makes the result depend on row order, and can flip which
# neighbouring bin gets flagged as altered.
# Bins uncovered by the truth default to diploid (1).
truth_cn_vec <- function(truth_dt, grid, node_name) {
    n_bins <- nrow(grid)
    v <- integer(2L * n_bins) + 1L
    sub <- truth_dt[name == node_name]
    if (nrow(sub) == 0L) return(NULL)

    setkey(grid, chrom, start, end)
    sub <- copy(sub)
    setnames(sub, c("start", "stop"), c("t_start", "t_end"))
    setkey(sub, chrom, t_start, t_end)
    ov <- foverlaps(sub, grid, by.x = c("chrom", "t_start", "t_end"),
                    type = "any", nomatch = 0L)

    # Overlap length in bp (approximate: half-open MEDICC2 vs closed truth,
    # 1-bp accuracy is fine here). max(0, .) not needed because foverlaps
    # already dropped non-overlapping rows.
    ov[, overlap := pmin(end, t_end) - pmax(start, t_start)]
    # Keep the truth segment with the largest overlap for each (bin, hap).
    winners <- ov[ov[, .I[which.max(overlap)], by = .(bin_id, haplotype)]$V1]

    col_offset <- ifelse(winners$haplotype == 1L, 0L, n_bins)
    v[winners$bin_id + col_offset] <- winners$cn
    v
}

# --- Metrics ---------------------------------------------------------------
# `t`, `m` are integer vectors of length 2 * n_bins. Returns a one-row tibble.
mrca_metrics <- function(t, m) {
    dip <- 1L
    stopifnot(length(t) == length(m))

    n_hbins <- length(t)
    t_alt   <- t != dip
    m_alt   <- m != dip
    disagree <- t != m

    truth_altered_hbins <- sum(t_alt)
    med_altered_hbins   <- sum(m_alt)
    truth_altered_units <- sum(abs(t - dip))
    med_altered_units   <- sum(abs(m - dip))

    disagree_hbins <- sum(disagree)
    disagree_units <- sum(abs(t - m))

    tp_hbins <- sum(t_alt & m_alt)          # both altered
    fp_hbins <- sum(!t_alt & m_alt)         # hallucination
    fn_hbins <- sum(t_alt & !m_alt)         # missed
    tn_hbins <- sum(!t_alt & !m_alt)        # both diploid

    # Hallucination severity: sum of MEDICC2 deviation from diploid at
    # truth-diploid hbins.
    halluc_mask <- !t_alt & m_alt
    halluc_units       <- sum(abs(m[halluc_mask] - dip))
    halluc_gain_hbins  <- sum(halluc_mask & (m > dip))
    halluc_loss_hbins  <- sum(halluc_mask & (m < dip))
    halluc_gain_units  <- sum((m[halluc_mask & (m > dip)] - dip))
    halluc_loss_units  <- sum((dip - m[halluc_mask & (m < dip)]))

    missed_mask  <- t_alt & !m_alt
    missed_units <- sum(abs(t[missed_mask] - dip))

    miscall_mask  <- t_alt & m_alt & disagree
    miscall_hbins <- sum(miscall_mask)
    miscall_units <- sum(abs(t[miscall_mask] - m[miscall_mask]))

    tibble(
        n_hbins = n_hbins,
        truth_is_diploid    = (truth_altered_hbins == 0L),
        med_is_diploid      = (med_altered_hbins == 0L),
        exact_match         = (disagree_hbins == 0L),
        truth_altered_hbins = truth_altered_hbins,
        truth_altered_units = truth_altered_units,
        med_altered_hbins   = med_altered_hbins,
        med_altered_units   = med_altered_units,
        disagree_hbins      = disagree_hbins,
        disagree_units      = disagree_units,
        tp_hbins = tp_hbins, fp_hbins = fp_hbins,
        fn_hbins = fn_hbins, tn_hbins = tn_hbins,
        halluc_hbins      = fp_hbins,
        halluc_units      = halluc_units,
        halluc_gain_hbins = halluc_gain_hbins,
        halluc_loss_hbins = halluc_loss_hbins,
        halluc_gain_units = halluc_gain_units,
        halluc_loss_units = halluc_loss_units,
        missed_hbins      = fn_hbins,
        missed_units      = missed_units,
        miscall_hbins     = miscall_hbins,
        miscall_units     = miscall_units
    )
}

# --- Per-sim driver --------------------------------------------------------
process_sim <- function(sim_rel) {
    truth_dir <- file.path(ROOT, "data", "CN_subsampled", dirname(sim_rel))
    med_dir   <- file.path(ROOT, "data", "MEDICC2", "treeinference", sim_rel)

    # Trees: collapse unary chains in the truth tree so its root corresponds
    # to the MRCA of the sampled cells (matching what MEDICC2 reconstructs).
    truth_raw <- read_truth_tree(truth_dir)
    med_full  <- read_medicc2_tree(med_dir)
    med_tree  <- if ("diploid" %in% med_full$tip.label)
        ape::drop.tip(med_full, "diploid") else med_full

    collapsed  <- collapse_unary(truth_raw)
    truth_tree <- collapsed$tree

    # MRCA labels: index 1 of internal-node labels is the root under
    # ape::prop.part ordering.
    truth_mrca_name <- truth_tree$node.label[1L]
    med_mrca_name   <- med_tree$node.label[1L]

    # Profiles
    med_dt   <- fread(file.path(med_dir,   "sim1_final_cn_profiles.tsv"))
    truth_dt <- fread(file.path(truth_dir, "sim1_truth_profiles.tsv"))

    grid <- build_bin_grid(med_dt)

    m_vec <- medicc_cn_vec(med_dt, grid, med_mrca_name)
    t_vec <- truth_cn_vec(truth_dt, grid, truth_mrca_name)
    if (is.null(m_vec) || is.null(t_vec)) {
        message(sprintf("  [%s] missing profile for MRCA -- skipping", sim_rel))
        return(NULL)
    }

    stats <- mrca_metrics(t_vec, m_vec)

    bind_cols(
        parse_sim_params(sim_rel),
        tibble(truth_mrca_name = truth_mrca_name,
               medicc_mrca_name = med_mrca_name),
        stats
    )
}

# --- Main loop -------------------------------------------------------------
existing <- if (file.exists(OUT_CSV) && !force) {
    suppressMessages(read_csv(OUT_CSV, show_col_types = FALSE))
} else tibble()

todo <- if (force) sim_rels else {
    have <- if (nrow(existing) > 0L) existing$sim_rel else character()
    sort(setdiff(sim_rels, have))
}
message(sprintf("todo: %d sims (of %d discovered)", length(todo), length(sim_rels)))

run_one <- function(sim_rel) {
    t0 <- Sys.time()
    res <- tryCatch(process_sim(sim_rel),
                    error = function(e) {
                        message(sprintf("  FAILED %s: %s", sim_rel, e$message))
                        NULL
                    })
    dt <- as.numeric(Sys.time() - t0, units = "secs")
    message(sprintf("  done %s (%.1fs)", sim_rel, dt))
    res
}

t0 <- Sys.time()
new_rows <- if (cores > 1L) {
    mclapply(todo, run_one, mc.cores = cores)
} else {
    lapply(todo, run_one)
}
dt <- as.numeric(Sys.time() - t0, units = "secs")
message(sprintf("processed %d sims in %.1fs (avg %.2fs/sim)",
                length(todo), dt, dt / max(length(todo), 1L)))

new_df <- bind_rows(Filter(Negate(is.null), new_rows))
combined <- if (force || nrow(existing) == 0L) {
    new_df
} else {
    kept <- existing %>% filter(!(sim_rel %in% new_df$sim_rel))
    bind_rows(kept, new_df)
}
combined <- combined %>% arrange(scenario, timing, N, d, s, M_inj)

tmp <- paste0(OUT_CSV, ".tmp")
write_csv(combined, tmp)
file.rename(tmp, OUT_CSV)
message(sprintf("wrote: %s (%d rows)", OUT_CSV, nrow(combined)))
