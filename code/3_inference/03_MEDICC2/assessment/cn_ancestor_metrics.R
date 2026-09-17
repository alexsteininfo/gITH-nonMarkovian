# cn_ancestor_metrics.R
#
# For each MEDICC2 sim on disk, compare the CN profiles MEDICC2 reconstructs at
# internal (ancestral) nodes against the truth CN profiles carried in the
# subsampled truth tree. Emits two tiers of output:
#
#   1. Per-sim per-node CSVs under data/MEDICC2/benchmark/cn_ancestor/<slug>.csv
#      -- one row per matched-clade (strategy A), per leaf-pair MRCA
#      (strategy B), and one row for the sample-MRCA. Six metrics per row:
#      Hamming and Sum-|Delta|, each raw / signal-restricted / baseline-relative.
#
#   2. A summary CSV data/MEDICC2/benchmark/cn_ancestor_summary.csv
#      -- one row per sim, with sample-MRCA values + per-strategy medians +
#      strategy-A matched fraction. Aligned to aggregate_metrics.csv keys.
#
# Two matching strategies (see brainstorm notes):
#   A. Clade-based (bipartition). Match a MEDICC2 internal node to a truth
#      internal node iff both subtend the same set of sampled leaves. Only
#      compares topologically-agreed nodes -- honest but selection-biased.
#   B. Pair-MRCA. For every unordered pair (i, j) of sampled leaves, take the
#      MRCA in each tree and compare its CN profile. Always defined; every
#      pair contributes one row.
#
# Six metrics per node comparison:
#   hamming_raw       -- fraction of segments where truth != MEDICC2
#   hamming_signal    -- disagreeing altered segments / union of altered segs
#                        (altered := CN != 2 on that haplotype)
#   hamming_baseline  -- 1 - sum(truth!=med) / sum(truth!=diploid)
#                        R^2-like against the "always diploid" baseline
#   sumabs_raw        -- sum |truth - MEDICC2| over segments (both haplotypes)
#   sumabs_signal     -- sum |truth - med| / sum max(|truth-2|, |med-2|)
#   sumabs_baseline   -- 1 - sum(|truth - med|) / sum(|truth - 2|)
#
# CN grid: MEDICC2 uses fixed 1 Mb bins; truth uses variable-length segments.
# We project truth onto MEDICC2's bin grid once per sim (interval overlap),
# yielding a length-(2 * n_bins) integer vector per node with haplotype A bins
# followed by haplotype B bins.
#
# Idempotent: per-sim CSVs already on disk are skipped unless --force.
# Optionally parallelised across sims via --cores N (default 1).
#
# Usage:
#   Rscript code/3_inference/03_MEDICC2/assessment/cn_ancestor_metrics.R \
#       [--force] [--cores N] [--sim <sim_rel>]

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(readr); library(tibble)
    library(ape); library(data.table); library(parallel)
})

args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE))
if (length(this_file) != 1L) stop("cannot resolve --file")
ROOT <- normalizePath(file.path(dirname(this_file), "..", "..", "..", ".."))

source(file.path(ROOT, "code", "X_helpers", "medicc2_io.R"))

# -------------------------- Arg parsing --------------------------------------
user_args <- commandArgs(trailingOnly = TRUE)
force <- FALSE
cores <- 1L
sim_filter <- NULL
i <- 1L
while (i <= length(user_args)) {
    a <- user_args[i]
    if (a == "--force") { force <- TRUE; i <- i + 1L }
    else if (a == "--cores") { cores <- as.integer(user_args[i + 1L]); i <- i + 2L }
    else if (a == "--sim")   { sim_filter <- user_args[i + 1L]; i <- i + 2L }
    else stop("unknown arg: ", a)
}
message(sprintf("force=%s cores=%d sim_filter=%s",
                force, cores, if (is.null(sim_filter)) "<all>" else sim_filter))

OUT_DIR    <- file.path(ROOT, "data", "MEDICC2", "benchmark", "cn_ancestor")
SUMMARY_CSV <- file.path(ROOT, "data", "MEDICC2", "benchmark", "cn_ancestor_summary.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MED_ROOT <- file.path(ROOT, "data", "MEDICC2", "treeinference")

# -------------------------- Sim discovery ------------------------------------
tree_files <- list.files(MED_ROOT, pattern = "^sim.*_final_tree\\.new$",
                         recursive = TRUE, full.names = TRUE)
sim_rels <- sort(unique(dirname(sub(paste0("^", MED_ROOT, "/"), "", tree_files))))
if (!is.null(sim_filter)) sim_rels <- intersect(sim_rels, sim_filter)
message(sprintf("Discovered %d sims", length(sim_rels)))

# Parameter parser (same shape as aggregate_metrics.R)
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

sim_csv_path <- function(sim_rel) {
    file.path(OUT_DIR, paste0(sim_slug(sim_rel), ".csv"))
}

# -------------------------- CN matrix builders -------------------------------
# Build a bin-grid data.frame (chrom, start, end, bin_id) from a MEDICC2
# profiles data.table -- the grid is the same across all samples in one sim.
build_bin_grid <- function(med_dt) {
    grid <- unique(med_dt[, .(chrom, start, end)])
    setorder(grid, chrom, start)
    grid[, bin_id := seq_len(.N)]
    grid[]
}

# MEDICC2 CN matrix: nodes as rows, 2*n_bins as columns (cn_a bins then cn_b bins).
build_medicc_cn_matrix <- function(med_dt, grid) {
    n_bins <- nrow(grid)
    med_dt <- merge(med_dt, grid, by = c("chrom", "start", "end"), sort = FALSE)
    nodes  <- sort(unique(med_dt$sample_id))
    node_idx <- setNames(seq_along(nodes), nodes)
    m <- matrix(NA_integer_, nrow = length(nodes), ncol = 2L * n_bins,
                dimnames = list(nodes, NULL))
    m[cbind(node_idx[med_dt$sample_id], med_dt$bin_id)]          <- med_dt$cn_a
    m[cbind(node_idx[med_dt$sample_id], med_dt$bin_id + n_bins)] <- med_dt$cn_b
    m
}

# Truth CN matrix, projected onto the MEDICC2 bin grid via interval overlap.
# haplotype 1 -> cn_a slot (columns 1..n_bins); haplotype 2 -> cn_b slot
# (columns n_bins+1..2*n_bins).
build_truth_cn_matrix <- function(truth_dt, grid) {
    n_bins <- nrow(grid)
    # data.table foverlaps: treat truth segments as [start, stop] and MEDICC2
    # bins as [start, end]. Overlap when they share at least 1 bp.
    setkey(grid, chrom, start, end)
    truth_dt <- copy(truth_dt)  # avoid mutating caller's dt
    setnames(truth_dt, c("start", "stop"), c("t_start", "t_end"))
    setkey(truth_dt, chrom, t_start, t_end)
    ov <- foverlaps(truth_dt, grid, by.x = c("chrom", "t_start", "t_end"),
                    type = "any", nomatch = 0L)
    # ov has one row per (truth segment, overlapping bin). For each such row
    # we assign truth's cn to the MEDICC2 bin.
    nodes <- sort(unique(truth_dt$name))
    node_idx <- setNames(seq_along(nodes), nodes)
    m <- matrix(2L, nrow = length(nodes), ncol = 2L * n_bins,
                dimnames = list(nodes, NULL))
    # Column offset by haplotype (1 -> a, 2 -> b)
    col_offset <- ifelse(ov$haplotype == 1L, 0L, n_bins)
    m[cbind(node_idx[ov$name], ov$bin_id + col_offset)] <- ov$cn
    m
}

# -------------------------- Metric computation -------------------------------
# `t` and `m` are integer vectors of length 2 * n_bins. Returns a named list.
cn_metrics <- function(t, m) {
    dip <- 2L
    diff  <- t != m
    n_seg <- length(t)

    hamming_raw <- sum(diff) / n_seg
    sumabs_raw  <- sum(abs(t - m))

    # Signal-restricted: at least one is altered
    altered <- (t != dip) | (m != dip)
    n_alt <- sum(altered)
    hamming_signal <- if (n_alt > 0L) sum(diff & altered) / n_alt else NA_real_

    denom_sumabs_sig <- sum(pmax(abs(t - dip), abs(m - dip)))
    sumabs_signal <- if (denom_sumabs_sig > 0L) sum(abs(t - m)) / denom_sumabs_sig else NA_real_

    # Baseline-relative (R^2-like against all-diploid prediction)
    denom_ham_base <- sum(t != dip)
    hamming_baseline <- if (denom_ham_base > 0L) 1 - sum(diff) / denom_ham_base else NA_real_

    denom_sumabs_base <- sum(abs(t - dip))
    sumabs_baseline <- if (denom_sumabs_base > 0L) 1 - sum(abs(t - m)) / denom_sumabs_base else NA_real_

    list(
        hamming_raw = hamming_raw, hamming_signal = hamming_signal,
        hamming_baseline = hamming_baseline,
        sumabs_raw = as.numeric(sumabs_raw), sumabs_signal = sumabs_signal,
        sumabs_baseline = sumabs_baseline
    )
}

metrics_row <- function(t_vec, m_vec) as_tibble(cn_metrics(t_vec, m_vec))

# -------------------------- Per-sim driver -----------------------------------
process_sim <- function(sim_rel) {
    truth_dir <- file.path(ROOT, "data", "CN_subsampled", dirname(sim_rel))
    med_dir   <- file.path(ROOT, "data", "MEDICC2", "treeinference", sim_rel)

    # Trees
    truth_raw <- read_truth_tree(truth_dir)
    med_full  <- read_medicc2_tree(med_dir)
    med_tree  <- if ("diploid" %in% med_full$tip.label)
        ape::drop.tip(med_full, "diploid") else med_full

    collapsed <- collapse_unary(truth_raw)
    truth_tree <- collapsed$tree

    common <- intersect(truth_tree$tip.label, med_tree$tip.label)
    if (length(common) < 0.9 * length(truth_tree$tip.label)) {
        message(sprintf("  [%s] <90%% tip overlap -- skipping", sim_rel))
        return(NULL)
    }
    # Restrict both trees to common tips (matches the RF/etc convention).
    if (length(common) < length(truth_tree$tip.label))
        truth_tree <- ape::drop.tip(truth_tree, setdiff(truth_tree$tip.label, common))
    if (length(common) < length(med_tree$tip.label))
        med_tree <- ape::drop.tip(med_tree, setdiff(med_tree$tip.label, common))

    # CN profiles -- data.table for the big MEDICC2 file
    med_dt <- fread(file.path(med_dir, "sim1_final_cn_profiles.tsv"))
    truth_dt <- fread(file.path(truth_dir, "sim1_truth_profiles.tsv"))

    grid   <- build_bin_grid(med_dt)
    n_bins <- nrow(grid)
    med_cn   <- build_medicc_cn_matrix(med_dt, grid)
    truth_cn <- build_truth_cn_matrix(truth_dt, grid)

    # Look up CN for a node label; NA if the label is not covered.
    get_cn <- function(mat, name) {
        idx <- match(name, rownames(mat))
        if (is.na(idx)) NULL else mat[idx, ]
    }

    # ----- Internal node labels + descendant leaf sets ----------------------
    truth_parts <- ape::prop.part(truth_tree)  # one entry per internal node, root first
    med_parts   <- ape::prop.part(med_tree)
    truth_int_labels <- truth_tree$node.label
    med_int_labels   <- med_tree$node.label
    stopifnot(length(truth_parts) == length(truth_int_labels))
    stopifnot(length(med_parts)   == length(med_int_labels))

    truth_leaf_sets <- lapply(truth_parts, function(idx) sort(truth_tree$tip.label[idx]))
    med_leaf_sets   <- lapply(med_parts,   function(idx) sort(med_tree$tip.label[idx]))
    hash <- function(x) paste(x, collapse = "|")
    truth_keys <- vapply(truth_leaf_sets, hash, character(1))
    med_keys   <- vapply(med_leaf_sets,   hash, character(1))
    clade_sizes_truth <- lengths(truth_leaf_sets)
    clade_sizes_med   <- lengths(med_leaf_sets)

    # ----- Strategy A: matched clades ---------------------------------------
    match_idx <- match(truth_keys, med_keys)
    matched_truth_pos <- which(!is.na(match_idx))
    matched_med_pos   <- match_idx[matched_truth_pos]

    rows_A <- lapply(seq_along(matched_truth_pos), function(k) {
        ti <- matched_truth_pos[k]; mi <- matched_med_pos[k]
        t_name <- truth_int_labels[ti]; m_name <- med_int_labels[mi]
        t_vec <- get_cn(truth_cn, t_name)
        m_vec <- get_cn(med_cn,   m_name)
        if (is.null(t_vec) || is.null(m_vec)) return(NULL)
        r <- metrics_row(t_vec, m_vec)
        tibble(
            strategy = "A", leaf_i = NA_character_, leaf_j = NA_character_,
            truth_mrca_name = t_name, medicc_mrca_name = m_name,
            clade_size_truth = clade_sizes_truth[ti],
            clade_size_medicc = clade_sizes_med[mi]
        ) %>% bind_cols(r)
    })
    df_A <- if (length(rows_A)) bind_rows(rows_A) else tibble()

    # ----- Strategy B: pair-MRCA --------------------------------------------
    n_tips <- length(common)
    # Precompute MRCA matrices (indexed by node position in each tree)
    truth_mrca_mat <- ape::mrca(truth_tree, full = FALSE)  # Ntip x Ntip
    med_mrca_mat   <- ape::mrca(med_tree,   full = FALSE)
    # tip label -> tip index in each tree
    truth_tip_idx <- setNames(seq_len(n_tips), truth_tree$tip.label)
    med_tip_idx   <- setNames(seq_len(n_tips), med_tree$tip.label)
    # Internal node index -> label + leaf-set-size lookups (ape indexes internal
    # nodes as (Ntip+1)..(2*Ntip-1) for a binary tree). prop.part order is the
    # same, so int_label_by_idx[k] gives the label of internal node Ntip+k.
    truth_int_by_node <- function(idx) truth_int_labels[idx - n_tips]
    med_int_by_node   <- function(idx) med_int_labels[idx - n_tips]
    truth_clade_by_node <- function(idx) clade_sizes_truth[idx - n_tips]
    med_clade_by_node   <- function(idx) clade_sizes_med[idx - n_tips]

    # Enumerate all pairs
    pair_mat <- t(combn(common, 2L))
    n_pairs <- nrow(pair_mat)
    ti_i <- truth_tip_idx[pair_mat[, 1L]]
    ti_j <- truth_tip_idx[pair_mat[, 2L]]
    mi_i <- med_tip_idx[pair_mat[, 1L]]
    mi_j <- med_tip_idx[pair_mat[, 2L]]
    t_mrca <- truth_mrca_mat[cbind(ti_i, ti_j)]
    m_mrca <- med_mrca_mat[cbind(mi_i, mi_j)]

    # Cache CN vectors per unique MRCA node so we don't recompute per-pair
    unique_t_mrca <- unique(t_mrca); unique_m_mrca <- unique(m_mrca)
    t_names_by_node <- setNames(truth_int_by_node(unique_t_mrca), unique_t_mrca)
    m_names_by_node <- setNames(med_int_by_node(unique_m_mrca),   unique_m_mrca)
    t_cache <- setNames(lapply(t_names_by_node, function(nm) get_cn(truth_cn, nm)),
                        unique_t_mrca)
    m_cache <- setNames(lapply(m_names_by_node, function(nm) get_cn(med_cn, nm)),
                        unique_m_mrca)

    rows_B <- lapply(seq_len(n_pairs), function(p) {
        t_vec <- t_cache[[as.character(t_mrca[p])]]
        m_vec <- m_cache[[as.character(m_mrca[p])]]
        if (is.null(t_vec) || is.null(m_vec)) return(NULL)
        r <- metrics_row(t_vec, m_vec)
        tibble(
            strategy = "B", leaf_i = pair_mat[p, 1L], leaf_j = pair_mat[p, 2L],
            truth_mrca_name = t_names_by_node[[as.character(t_mrca[p])]],
            medicc_mrca_name = m_names_by_node[[as.character(m_mrca[p])]],
            clade_size_truth = truth_clade_by_node(t_mrca[p]),
            clade_size_medicc = med_clade_by_node(m_mrca[p])
        ) %>% bind_cols(r)
    })
    df_B <- if (length(rows_B)) bind_rows(rows_B) else tibble()

    # ----- Sample-MRCA (root of each tree) ----------------------------------
    truth_root_label <- truth_int_labels[1L]
    med_root_label   <- med_int_labels[1L]
    smr_row <- {
        t_vec <- get_cn(truth_cn, truth_root_label)
        m_vec <- get_cn(med_cn,   med_root_label)
        if (is.null(t_vec) || is.null(m_vec)) tibble() else {
            r <- metrics_row(t_vec, m_vec)
            tibble(
                strategy = "sample_mrca",
                leaf_i = NA_character_, leaf_j = NA_character_,
                truth_mrca_name = truth_root_label,
                medicc_mrca_name = med_root_label,
                clade_size_truth = n_tips, clade_size_medicc = n_tips
            ) %>% bind_cols(r)
        }
    }

    combined <- bind_rows(df_A, df_B, smr_row)

    # ----- Persist per-sim CSV atomically -----------------------------------
    out_path <- sim_csv_path(sim_rel)
    tmp_path <- paste0(out_path, ".tmp")
    fwrite(combined, tmp_path)
    file.rename(tmp_path, out_path)

    # ----- Summary row ------------------------------------------------------
    med_A <- if (nrow(df_A) > 0L) df_A else tibble()
    med_B <- if (nrow(df_B) > 0L) df_B else tibble()
    med_or_NA <- function(x) if (length(x) > 0L) median(x, na.rm = TRUE) else NA_real_
    smr <- if (nrow(smr_row) > 0L) smr_row else NULL

    sample_val <- function(col) if (!is.null(smr)) smr[[col]] else NA_real_
    med_val    <- function(df, col) if (nrow(df) > 0L) med_or_NA(df[[col]]) else NA_real_

    tibble(
        A_n_matched      = nrow(df_A),
        A_n_medicc_total = length(med_int_labels),
        B_n_pairs        = nrow(df_B),
        # sample-MRCA (headline)
        sample_mrca_hamming_raw      = sample_val("hamming_raw"),
        sample_mrca_hamming_signal   = sample_val("hamming_signal"),
        sample_mrca_hamming_baseline = sample_val("hamming_baseline"),
        sample_mrca_sumabs_raw       = sample_val("sumabs_raw"),
        sample_mrca_sumabs_signal    = sample_val("sumabs_signal"),
        sample_mrca_sumabs_baseline  = sample_val("sumabs_baseline"),
        # Strategy A medians
        A_med_hamming_raw      = med_val(df_A, "hamming_raw"),
        A_med_hamming_signal   = med_val(df_A, "hamming_signal"),
        A_med_hamming_baseline = med_val(df_A, "hamming_baseline"),
        A_med_sumabs_raw       = med_val(df_A, "sumabs_raw"),
        A_med_sumabs_signal    = med_val(df_A, "sumabs_signal"),
        A_med_sumabs_baseline  = med_val(df_A, "sumabs_baseline"),
        # Strategy B medians
        B_med_hamming_raw      = med_val(df_B, "hamming_raw"),
        B_med_hamming_signal   = med_val(df_B, "hamming_signal"),
        B_med_hamming_baseline = med_val(df_B, "hamming_baseline"),
        B_med_sumabs_raw       = med_val(df_B, "sumabs_raw"),
        B_med_sumabs_signal    = med_val(df_B, "sumabs_signal"),
        B_med_sumabs_baseline  = med_val(df_B, "sumabs_baseline")
    )
}

# -------------------------- Main loop ----------------------------------------
existing_summary <- if (file.exists(SUMMARY_CSV) && !force) {
    suppressMessages(read_csv(SUMMARY_CSV, show_col_types = FALSE))
} else tibble()

todo <- if (force) sim_rels else {
    have <- if (nrow(existing_summary) > 0L) existing_summary$sim_rel else character()
    also_missing_csv <- setdiff(sim_rels,
                                 sim_rels[file.exists(vapply(sim_rels, sim_csv_path, character(1)))])
    sort(unique(c(setdiff(sim_rels, have), also_missing_csv)))
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
    if (is.null(res)) NULL else bind_cols(parse_sim_params(sim_rel), res)
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

combined <- if (force || nrow(existing_summary) == 0L) {
    new_df
} else {
    kept <- existing_summary %>% filter(!(sim_rel %in% new_df$sim_rel))
    bind_rows(kept, new_df)
}
combined <- combined %>% arrange(scenario, timing, N, d, s, M_inj)

tmp <- paste0(SUMMARY_CSV, ".tmp")
write_csv(combined, tmp)
file.rename(tmp, SUMMARY_CSV)
message(sprintf("wrote: %s (%d rows)", SUMMARY_CSV, nrow(combined)))
