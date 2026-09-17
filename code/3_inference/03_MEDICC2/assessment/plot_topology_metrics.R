# plot_topology_metrics.R
#
# Tree topology + branch-length agreement between the ground-truth
# (subsampled, unary-collapsed) tree and the MEDICC2 reconstruction, for one
# simulation. Implements the four metrics recommended in
# `tree_similarity_metrics.md`:
#   1. Normalised Robinson-Foulds (topology only)
#   2. Weighted Robinson-Foulds (topology + branch length in events)
#   3. Cophenetic correlation (pairwise-distance rank + linear slope)
#   4. Sister-pair (cherry) precision / recall / F1 -- the pure-leaf analogue
#      of ancestor-descendant F1 (the AD matrix collapses to trivial on a
#      leaves-only tree; sister pairs are the equivalent local-relationship
#      signal).
#
# Truth-tree branch lengths are assigned as *events per collapsed edge*
# (from `edge_event_counts_truth`) so that both trees share the "events"
# unit; this is the only assignment that makes weighted RF and cophenetic
# distances directly comparable to MEDICC2's output.

suppressPackageStartupMessages({
    library(dplyr); library(readr); library(ggplot2)
    library(ape); library(phangorn); library(cowplot); library(gridExtra); library(grid)
})

# ROOT resolution + arg parsing (matches the other plot_*.R scripts)
args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE))
if (length(this_file) != 1L) stop("cannot resolve --file")
ROOT <- normalizePath(file.path(dirname(this_file), "..", "..", "..", ".."))

source(file.path(ROOT, "code", "X_helpers", "medicc2_io.R"))
source(file.path(ROOT, "code", "X_helpers", "plotting_functions.R"))

user_args <- commandArgs(trailingOnly = TRUE)
sim_rel <- "neutral/gamma/neutral_gamma_N1000_d0.5_k5.0_n100/sim1"
outdir_arg <- NULL
i <- 1
while (i <= length(user_args)) {
    if (user_args[i] == "--sim") { sim_rel <- user_args[i + 1L]; i <- i + 2L }
    else if (user_args[i] == "--outdir") { outdir_arg <- user_args[i + 1L]; i <- i + 2L }
    else stop("unknown arg: ", user_args[i])
}
message("sim: ", sim_rel)

TRUTH_DIR   <- file.path(ROOT, "data", "CN_subsampled", dirname(sim_rel))
MEDICC2_DIR <- file.path(ROOT, "data", "MEDICC2", "treeinference", sim_rel)
OUT_DIR <- if (is.null(outdir_arg)) {
    file.path(ROOT, "figures", "3_inference", "03_MEDICC2", "parsimony_saturation", "cophenetic")
} else if (startsWith(outdir_arg, "/")) {
    outdir_arg
} else {
    file.path(ROOT, outdir_arg)
}
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Load and prepare trees so both use the "events" branch-length unit ---
truth_tree_raw    <- read_truth_tree(TRUTH_DIR)
truth_events      <- read_truth_events(TRUTH_DIR)
medicc2_tree_full <- read_medicc2_tree(MEDICC2_DIR)

collapsed  <- collapse_unary(truth_tree_raw)
truth_tree <- collapsed$tree
ec         <- edge_event_counts_truth(collapsed, truth_events)
# ec$edge_id runs 1..nrow(truth_tree$edge) in tree$edge order per the helper contract
truth_tree$edge.length <- ec$events

if ("diploid" %in% medicc2_tree_full$tip.label) {
    medicc2_tree <- ape::drop.tip(medicc2_tree_full, "diploid")
} else {
    medicc2_tree <- medicc2_tree_full
}

common <- intersect(truth_tree$tip.label, medicc2_tree$tip.label)
message(sprintf("common tips: %d (truth=%d, medicc2=%d)",
                length(common), length(truth_tree$tip.label),
                length(medicc2_tree$tip.label)))
if (length(common) < 0.9 * length(truth_tree$tip.label)) {
    stop("<90% tip overlap between truth and MEDICC2 -- refusing to compute metrics")
}
truth_c <- if (length(common) < length(truth_tree$tip.label))
    ape::drop.tip(truth_tree, setdiff(truth_tree$tip.label, common)) else truth_tree
med_c   <- if (length(common) < length(medicc2_tree$tip.label))
    ape::drop.tip(medicc2_tree, setdiff(medicc2_tree$tip.label, common)) else medicc2_tree

# --- 1 & 2. Robinson-Foulds (topology only, and length-weighted in events) ---
rf_raw   <- phangorn::RF.dist(truth_c, med_c, normalize = FALSE)
rf_norm  <- phangorn::RF.dist(truth_c, med_c, normalize = TRUE)
wrf_raw  <- phangorn::wRF.dist(truth_c, med_c, normalize = FALSE)
wrf_norm <- phangorn::wRF.dist(truth_c, med_c, normalize = TRUE)

# --- 3. Cophenetic pairwise distances (events on the path between cell pairs) ---
co_t <- ape::cophenetic.phylo(truth_c)[common, common]
co_m <- ape::cophenetic.phylo(med_c)[common, common]
uti  <- upper.tri(co_t)
pair_df <- tibble(truth = co_t[uti], medicc2 = co_m[uti])

r_pearson  <- cor(pair_df$truth, pair_df$medicc2, method = "pearson")
r_spearman <- cor(pair_df$truth, pair_df$medicc2, method = "spearman")
fit <- lm(medicc2 ~ truth, data = pair_df)
slope     <- unname(coef(fit)["truth"])
intercept <- unname(coef(fit)["(Intercept)"])

# --- 4. Sister-pair (cherry) F1 ---
sister_pairs <- function(tr) {
    tips <- tr$tip.label
    n <- length(tips)
    children_of <- split(tr$edge[, 2], tr$edge[, 1])
    pairs <- lapply(children_of, function(kids) {
        # A cherry is an internal node with exactly two children that are both leaves.
        if (length(kids) == 2L && all(kids <= n)) sort(tips[kids]) else NULL
    })
    pairs <- Filter(Negate(is.null), pairs)
    if (!length(pairs)) return(character(0))
    vapply(pairs, function(p) paste(p, collapse = "|"), character(1))
}
s_truth <- sister_pairs(truth_c)
s_med   <- sister_pairs(med_c)
tp <- length(intersect(s_truth, s_med))
fp <- length(setdiff(s_med, s_truth))
fn <- length(setdiff(s_truth, s_med))
sister_prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
sister_rec  <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
sister_f1   <- if (is.finite(sister_prec + sister_rec) && (sister_prec + sister_rec) > 0)
    2 * sister_prec * sister_rec / (sister_prec + sister_rec) else NA_real_

# --- Panels ---
col_truth <- "#6BAED6"; col_med <- "#08306B"

p_scatter <- ggplot(pair_df, aes(truth, medicc2)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dashed") +
    geom_point(alpha = 0.2, size = 0.6, colour = col_med) +
    geom_smooth(method = "lm", colour = "#FB6A4A", fill = "#FCBBA1",
                se = TRUE, formula = y ~ x, linewidth = 0.6) +
    coord_equal() +
    labs(title = "Cophenetic pairwise distances (events on the path between cell pairs)",
         subtitle = sprintf(
             "Pearson r=%.3f, Spearman r=%.3f  |  slope=%.3f, intercept=%.2f  |  dashed = y=x",
             r_pearson, r_spearman, slope, intercept),
         x = "Truth: events on the path between the two cells",
         y = "MEDICC2: events on the path between the two cells") +
    theme_bw(base_size = 11)

fmt <- function(x) formatC(x, digits = 3, format = "g")
metrics_tbl <- tibble(
    Metric = c(
        "Robinson-Foulds (raw)",
        "Robinson-Foulds (normalised, 0..1)",
        "Weighted RF (raw, events)",
        "Weighted RF (normalised)",
        "Cophenetic r (Pearson)",
        "Cophenetic r (Spearman)",
        "Cophenetic slope (MEDICC2 ~ truth)",
        "Sister-pair precision",
        "Sister-pair recall",
        "Sister-pair F1"
    ),
    Value = c(
        fmt(rf_raw), fmt(rf_norm),
        fmt(wrf_raw), fmt(wrf_norm),
        fmt(r_pearson), fmt(r_spearman),
        fmt(slope),
        fmt(sister_prec), fmt(sister_rec), fmt(sister_f1)
    ),
    Interpretation = c(
        "splits present in one tree but not the other; lower = more agreement",
        "raw / (2n - 6); 0 = identical topology, 1 = maximally different",
        "sum of |edge-length differences| over matched splits, in events",
        "raw / max possible; smaller = better length agreement",
        "correlation of pairwise event distances; 1 = perfect",
        "rank correlation; 1 = ranking preserved",
        "slope of MEDICC2 pair-distances vs truth; <1 = branches compressed",
        "TP / (TP + FP) on cherries (sister-leaf pairs)",
        "TP / (TP + FN) on cherries",
        "harmonic mean of sister precision and recall"
    )
)

p_table <- gridExtra::tableGrob(
    metrics_tbl, rows = NULL,
    theme = gridExtra::ttheme_minimal(base_size = 9, padding = unit(c(4, 3), "mm"))
)

title_grob <- cowplot::ggdraw() + cowplot::draw_label(
    paste0("Tree topology + branch-length metrics\n", sim_rel),
    fontface = "bold", size = 12
)

fig <- cowplot::plot_grid(
    title_grob, p_scatter, p_table,
    ncol = 1, rel_heights = c(0.06, 0.52, 0.42)
)

out_path <- file.path(OUT_DIR, paste0("topology_metrics__", sim_slug(sim_rel), ".png"))
ggsave(out_path, fig, width = 9, height = 12, dpi = 150)
message("wrote: ", out_path)
message("")
message("Metrics:")
print(metrics_tbl)
