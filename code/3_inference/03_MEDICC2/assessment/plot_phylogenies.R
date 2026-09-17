# plot_phylogenies.R
#
# (B) Truth-vs-MEDICC2 phylogeny comparison. Writes:
#   phylogenies_paired__<slug>.png     -- two ggtree panels, matched tip order
#   phylogenies_tanglegram__<slug>.png -- ape::cophyloplot connecting lines
#
# Usage:
#   Rscript code/3_inference/03_MEDICC2/assessment/plot_phylogenies.R \
#       [--sim <sim_rel>] [--outdir <abs_path_or_repo_rel>]
#
# Inputs:
#   data/CN_subsampled/<dirname(sim_rel)>/sim1_tree.nwk
#   data/MEDICC2/treeinference/<sim_rel>/sim1_final_tree.new
#
# Outputs (default outdir = figures/3_inference/03_MEDICC2/phylogenies/):
#   <outdir>/paired/phylogenies_paired__<slug>.png
#   <outdir>/tanglegram/phylogenies_tanglegram__<slug>.png

suppressPackageStartupMessages({
    library(dplyr); library(ape); library(ggtree); library(ggplot2); library(cowplot)
})

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
# Absolute --outdir wins; a relative --outdir is resolved against ROOT; default
# is phylogenies/ (moved out of assessment/ to keep the tree comparisons grouped).
OUT_DIR <- if (is.null(outdir_arg)) {
    file.path(ROOT, "figures", "3_inference", "03_MEDICC2", "phylogenies")
} else if (startsWith(outdir_arg, "/")) {
    outdir_arg
} else {
    file.path(ROOT, outdir_arg)
}
dir.create(file.path(OUT_DIR, "paired"),    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "tanglegram"), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Step 2: Load trees, collapse truth, drop `diploid` from MEDICC2 for the visual
# ---------------------------------------------------------------------------

truth_tree_raw <- read_truth_tree(TRUTH_DIR)
medicc2_tree_full <- read_medicc2_tree(MEDICC2_DIR)

collapsed <- collapse_unary(truth_tree_raw)
truth_tree <- collapsed$tree

# Drop the diploid outgroup from the MEDICC2 tree for the visual only.
if ("diploid" %in% medicc2_tree_full$tip.label) {
    medicc2_tree <- ape::drop.tip(medicc2_tree_full, "diploid")
} else {
    medicc2_tree <- medicc2_tree_full
}

# Confirm shared tip set (should be the 100 sampled cells).
truth_tips <- truth_tree$tip.label
med_tips   <- medicc2_tree$tip.label
common <- intersect(truth_tips, med_tips)
message("Truth tips: ", length(truth_tips),
        "  MEDICC2 tips (no diploid): ", length(med_tips),
        "  Common: ", length(common))
if (length(common) < 0.9 * length(truth_tips)) {
    stop("plot_phylogenies: <90% tip overlap (truth=", length(truth_tips),
         ", medicc2=", length(med_tips), ", common=", length(common), ")")
}

# ---------------------------------------------------------------------------
# Step 3: Paired ggtree panels with matched leaf order
# ---------------------------------------------------------------------------

# Rotate the MEDICC2 tree to match the truth tip order as best we can.
# ape::rotateConstr tries to match tips in the given order.
med_tree_rot <- ape::rotateConstr(medicc2_tree, truth_tips)

# Colour tips by their truth-tree y-position (top -> bottom, viridis).
# Same tip label gets the same colour in both panels, so a smooth colour
# gradient on the MEDICC2 panel <=> tips are perfectly co-ordered with truth.
# Because rotateConstr can only rotate around internal nodes (never override
# topology), any disagreement between the two topologies shows up here as a
# region of "scrambled" colours in the MEDICC2 panel.
p_truth_base <- ggtree(truth_tree)
truth_tip_pos <- p_truth_base$data %>%
    dplyr::filter(isTip) %>%
    dplyr::transmute(label = label, tip_order = y)

viridis_scale <- ggplot2::scale_color_viridis_c(
    name = "Truth tip order (top -> bottom)",
    option = "viridis",
    guide = ggplot2::guide_colorbar(barwidth = 15, barheight = 0.6)
)

p_truth <- p_truth_base %<+% truth_tip_pos +
    ggtree::geom_tippoint(aes(color = tip_order), size = 1.4, na.rm = TRUE) +
    viridis_scale +
    ggtitle("Truth (collapsed lineage tree)") +
    theme_tree() +
    theme(legend.position = "none")

p_med <- ggtree(med_tree_rot) %<+% truth_tip_pos +
    ggtree::geom_tippoint(aes(color = tip_order), size = 1.4, na.rm = TRUE) +
    viridis_scale +
    ggtitle("MEDICC2 inferred") +
    theme_tree() +
    theme(legend.position = "none")

# Build a shared bottom legend so the colour bar isn't duplicated per panel.
legend_plot <- p_truth + theme(legend.position = "bottom")
shared_legend <- cowplot::get_legend(legend_plot)

paired_panels <- cowplot::plot_grid(p_truth, p_med, ncol = 2, align = "h")
paired <- cowplot::plot_grid(paired_panels, shared_legend,
                             ncol = 1, rel_heights = c(1, 0.06))
out_paired <- file.path(OUT_DIR, "paired", paste0("phylogenies_paired__", sim_slug(sim_rel), ".png"))
ggsave(out_paired, paired, width = 12, height = 10, dpi = 120)
message("wrote: ", out_paired)

# ---------------------------------------------------------------------------
# Step 4: Tanglegram via ape::cophyloplot
# ---------------------------------------------------------------------------

# assoc: 2-column matrix of matched tip names (only common tips).
# Uses the full medicc2_tree (without diploid already dropped above for visual).
assoc <- cbind(common, common)

out_tangle <- file.path(OUT_DIR, "tanglegram", paste0("phylogenies_tanglegram__", sim_slug(sim_rel), ".png"))
png(out_tangle, width = 12 * 120, height = 10 * 120, res = 120)
op <- par(mar = c(2, 1, 3, 1))
# Pass the rotated MEDICC2 tree (same object the paired plot uses). Because
# rotateConstr already rearranged children to match truth's tip order as
# closely as topology permits, connecting lines that remain diagonal here
# correspond to genuine topological disagreements that no rotation can fix.
ape::cophyloplot(truth_tree, med_tree_rot, assoc = assoc,
                 use.edge.length = FALSE,
                 show.tip.label = FALSE,
                 gap = 2, space = 40,
                 col = adjustcolor("black", alpha.f = 0.3))
title(main = paste0("Tanglegram - truth (left) vs MEDICC2 (right)\n", sim_rel),
      cex.main = 1.0)
par(op)
dev.off()
message("wrote: ", out_tangle)
