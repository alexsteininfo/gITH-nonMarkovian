# medicc2_io.R
#
# Shared loaders and tree-walking helpers for the MEDICC2 assessment
# scripts under code/3_inference/03_MEDICC2/assessment/.
#
# All loaders take an absolute directory path (`truth_dir` for CN_subsampled
# inputs, `medicc2_dir` for MEDICC2 outputs). Callers do path resolution.
#
# No plotting here -- palettes and themes live in plotting_functions.R.

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
    library(ape)
})

# Turn a repo-relative sim path like
#   "neutral/gamma/neutral_gamma_N1000_d0.5_k5.0_n100/sim1"
# into a filename-safe slug:
#   "neutral_gamma_N1000_d0.5_k5.0_n100__sim1"
# We use the last two components (the param dir + the sim number).
sim_slug <- function(sim_relpath) {
    parts <- strsplit(sim_relpath, "/", fixed = TRUE)[[1]]
    if (length(parts) < 2L) stop("sim_slug: expected at least 2 path components, got: ", sim_relpath)
    paste(tail(parts, 2), collapse = "__")
}

read_truth_tree <- function(truth_dir) {
    path <- file.path(truth_dir, "sim1_tree.nwk")
    if (!file.exists(path)) stop("read_truth_tree: not found: ", path)
    ape::read.tree(path)
}

read_medicc2_tree <- function(medicc2_dir) {
    path <- file.path(medicc2_dir, "sim1_final_tree.new")
    if (!file.exists(path)) stop("read_medicc2_tree: not found: ", path)
    ape::read.tree(path)
}

# MEDICC2 only tracks chr1..chr22, but the simulation truth includes chrX on
# both haplotypes (~4-5% of all events). Any statistic that will be compared
# against MEDICC2 output has to drop chrX first, otherwise truth-side counts
# are systematically inflated. Default is to drop chrX; pass keep_chrX = TRUE
# when analysing the truth in its own right.
read_truth_events <- function(truth_dir, keep_chrX = FALSE) {
    path <- file.path(truth_dir, "sim1_truth_events.tsv")
    ev <- read_tsv(path, col_types = cols(
        node_id = "i", name = "c", order = "i", type = "c",
        chrom = "c", haplotype = "i", start = "d", stop = "d",
        delta = "i", scale = "c", mode = "c"
    ))
    if (!keep_chrX) ev <- ev %>% filter(chrom != "chrX")
    ev
}

read_medicc2_events <- function(medicc2_dir) {
    path <- file.path(medicc2_dir, "sim1_copynumber_events_df.tsv")
    read_tsv(path, col_types = cols(
        sample_id = "c", chrom = "c", start = "d", end = "d",
        allele = "c", type = "c", cn_child = "i"
    ))
}

# Same chrX story as read_truth_events: MEDICC2 profiles only cover
# chr1..chr22, so callers that project truth onto MEDICC2's grid should drop
# chrX from the truth. (Projections via `foverlaps` against the MEDICC2 bin
# grid drop chrX by accident, but a caller iterating over truth rows directly
# would be off by ~5% otherwise.) Pass keep_chrX = TRUE for truth-only work.
read_truth_profiles <- function(truth_dir, keep_chrX = FALSE) {
    path <- file.path(truth_dir, "sim1_truth_profiles.tsv")
    prof <- read_tsv(path, col_types = cols(
        node_id = "i", name = "c", chrom = "c",
        haplotype = "i", start = "d", stop = "d", cn = "i"
    ))
    if (!keep_chrX) prof <- prof %>% filter(chrom != "chrX")
    prof
}

read_medicc2_profiles <- function(medicc2_dir) {
    # Large file (39 MB for n=100). read_tsv is fast enough.
    path <- file.path(medicc2_dir, "sim1_final_cn_profiles.tsv")
    read_tsv(path, col_types = cols(
        sample_id = "c", chrom = "c", start = "d", end = "d",
        cn_a = "i", cn_b = "i",
        is_gain = "l", is_loss = "l", is_wgd = "l",
        is_normal = "l", is_clonal = "l"
    ))
}

read_medicc2_branch_lengths <- function(medicc2_dir) {
    # Headerless 2-column TSV: <node>\t<length>
    path <- file.path(medicc2_dir, "sim1_branch_lengths.tsv")
    read_tsv(path, col_names = c("node", "length"), col_types = "cd")
}

# Return (tree, chain_nodes, divisions_collapsed) where each collapsed edge in
# `tree` carries the list of original node names on the unary chain strictly
# between its surviving parent (exclusive) and its surviving child (inclusive).
collapse_unary <- function(phy) {
    stopifnot(inherits(phy, "phylo"))
    if (is.null(phy$node.label)) {
        stop("collapse_unary: expected internal node labels on the truth tree")
    }

    all_names <- c(phy$tip.label, phy$node.label)
    n_tip <- length(phy$tip.label)
    n_node <- phy$Nnode

    # child -> parent map on the original tree
    parent_of <- integer(n_tip + n_node)
    parent_of[phy$edge[, 2]] <- phy$edge[, 1]
    # root has parent 0
    root_idx <- setdiff(phy$edge[, 1], phy$edge[, 2])[1]

    # A node "survives" iff it is a tip OR it has != 1 children (multifurcation
    # in a bifurcating tree = 2 children). Root survives by definition.
    n_children <- tabulate(phy$edge[, 1], nbins = n_tip + n_node)
    survives <- (seq_len(n_tip + n_node) <= n_tip) | (n_children != 1L)
    survives[root_idx] <- TRUE

    # Collapse via ape, then re-derive the chain nodes.
    tr <- ape::collapse.singles(phy)

    # For each edge in `tr`, walk from child upward on the ORIGINAL tree until we
    # hit the surviving parent; the walked nodes (child inclusive, parent
    # exclusive) form the chain.
    tr_all_names <- c(tr$tip.label, tr$node.label)
    orig_index <- match(tr_all_names, all_names)
    if (anyNA(orig_index)) {
        stop("collapse_unary: collapsed tree contains a name not in the original tree")
    }

    chain_nodes <- vector("list", nrow(tr$edge))
    for (i in seq_len(nrow(tr$edge))) {
        parent_new <- tr$edge[i, 1]
        child_new  <- tr$edge[i, 2]
        parent_orig <- orig_index[parent_new]
        child_orig  <- orig_index[child_new]

        chain <- character()
        cur <- child_orig
        while (cur != parent_orig) {
            chain <- c(all_names[cur], chain)  # prepend so order is root->leaf
            cur <- parent_of[cur]
            if (cur == 0L) {
                stop("collapse_unary: walked past root without hitting parent (edge ", i, ")")
            }
        }
        chain_nodes[[i]] <- chain
    }

    list(
        tree = tr,
        chain_nodes = chain_nodes,
        divisions_collapsed = vapply(chain_nodes, length, integer(1))
    )
}

# Count events per collapsed edge on the truth tree.
edge_event_counts_truth <- function(collapsed_res, events) {
    tr <- collapsed_res$tree
    chain_nodes <- collapsed_res$chain_nodes
    stopifnot(is.list(chain_nodes), length(chain_nodes) == nrow(tr$edge))

    events_by_name <- events %>%
        count(name, name = "n_events") %>%
        tibble::deframe()

    tr_all_names <- c(tr$tip.label, tr$node.label)

    tibble(
        edge_id = seq_along(chain_nodes),
        child = tr_all_names[tr$edge[, 2]],
        divisions = collapsed_res$divisions_collapsed,
        events = vapply(chain_nodes, function(ch) {
            sum(events_by_name[ch], na.rm = TRUE)
        }, numeric(1)) %>% as.integer()
    )
}

# Sum truth events on the root->leaf path in the collapsed truth tree.
per_cell_burden_truth <- function(collapsed_res, events) {
    tr <- collapsed_res$tree
    tips <- tr$tip.label

    ec <- edge_event_counts_truth(collapsed_res, events)  # tibble w/ edge_id, child, events

    # Map each tip to its ancestor edge chain (root-to-tip).
    # ape's nodepath gives node indices from root to tip.
    root_idx <- setdiff(tr$edge[, 1], tr$edge[, 2])[1]

    burden <- vapply(seq_along(tips), function(i) {
        tip_idx <- i  # tips are numbered 1..Ntip
        path <- ape::nodepath(tr, from = root_idx, to = tip_idx)
        # Each consecutive (parent, child) pair is an edge in tr$edge.
        edge_ids <- vapply(seq_len(length(path) - 1L), function(k) {
            which(tr$edge[, 1] == path[k] & tr$edge[, 2] == path[k + 1L])
        }, integer(1))
        sum(ec$events[edge_ids])
    }, numeric(1)) %>% as.integer()

    tibble(cell = tips, M = burden)
}

# Sum MEDICC2 branch lengths on the root->leaf path.
per_cell_burden_medicc2 <- function(tree, bl_df) {
    stopifnot(inherits(tree, "phylo"))
    tips <- tree$tip.label
    root_idx <- setdiff(tree$edge[, 1], tree$edge[, 2])[1]

    # Prefer tree$edge.length if it matches bl_df (MEDICC2 writes both).
    edge_len <- tree$edge.length
    if (is.null(edge_len)) {
        # Fallback: match branch lengths from bl_df by child node name.
        all_names <- c(tree$tip.label, tree$node.label)
        child_names <- all_names[tree$edge[, 2]]
        matched <- bl_df$length[match(child_names, bl_df$node)]
        if (anyNA(matched)) {
            missing <- child_names[is.na(matched)]
            stop("per_cell_burden_medicc2: branch_lengths.tsv missing entries for: ",
                 paste(head(missing, 5), collapse = ", "))
        }
        edge_len <- matched
    }

    burden <- vapply(seq_along(tips), function(i) {
        path <- ape::nodepath(tree, from = root_idx, to = i)
        edge_ids <- vapply(seq_len(length(path) - 1L), function(k) {
            which(tree$edge[, 1] == path[k] & tree$edge[, 2] == path[k + 1L])
        }, integer(1))
        sum(edge_len[edge_ids])
    }, numeric(1))

    tibble(cell = tips, M = burden)
}

# Pick n_stratified cells (min, 33rd pctile, 66th pctile, max by M) plus
# n_random uniform-random extras. Reproducible via `seed`.
stratified_cell_sample <- function(burden_df, n_stratified = 4L, n_random = 2L, seed = 42L) {
    stopifnot(all(c("cell", "M") %in% names(burden_df)), n_stratified >= 2L)

    df <- burden_df %>% arrange(M)
    n <- nrow(df)
    if (n < n_stratified + n_random) {
        stop("stratified_cell_sample: too few cells (", n, ") for ",
             n_stratified, " + ", n_random)
    }

    idxs_strat <- unique(round(seq(1, n, length.out = n_stratified)))
    strat_cells <- df$cell[idxs_strat]

    set.seed(seed)
    pool <- setdiff(df$cell, strat_cells)
    rand_cells <- if (n_random > 0L) sample(pool, n_random) else character()

    c(strat_cells, rand_cells)
}
