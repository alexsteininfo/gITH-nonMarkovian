# plot_saturation.R
#
# (C) Parsimony-saturation diagnostics: branch-length + per-cell burden
# histograms (truth vs MEDICC2) plus a moments table.

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(readr)
    library(ape); library(ggplot2); library(cowplot); library(gridExtra); library(grid)
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
OUT_DIR <- if (is.null(outdir_arg)) {
    file.path(ROOT, "figures", "3_inference", "03_MEDICC2", "parsimony_saturation")
} else if (startsWith(outdir_arg, "/")) {
    outdir_arg
} else {
    file.path(ROOT, outdir_arg)
}
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

truth_tree_raw <- read_truth_tree(TRUTH_DIR)
truth_events   <- read_truth_events(TRUTH_DIR)
medicc2_tree   <- read_medicc2_tree(MEDICC2_DIR)
bl_df          <- read_medicc2_branch_lengths(MEDICC2_DIR)

collapsed <- collapse_unary(truth_tree_raw)

# Branch lengths on both trees
ec_truth   <- edge_event_counts_truth(collapsed, truth_events)
bl_medicc2 <- bl_df$length

branch_df <- bind_rows(
    tibble(source = "truth",   len = as.numeric(ec_truth$events)),
    tibble(source = "MEDICC2", len = as.numeric(bl_medicc2))
)

# Per-cell burden (drop diploid on the MEDICC2 side for a like-for-like set)
bt <- per_cell_burden_truth(collapsed, truth_events)
bm <- per_cell_burden_medicc2(medicc2_tree, bl_df) %>%
    filter(cell != "diploid")

burden_df <- bind_rows(
    tibble(source = "truth",   cell = bt$cell, M = as.numeric(bt$M)),
    tibble(source = "MEDICC2", cell = bm$cell, M = as.numeric(bm$M))
)

# Moments
moment_row <- function(M) tibble(EM = mean(M), VarM = var(M),
                                 D = var(M) / mean(M),
                                 D_minus_1 = var(M) / mean(M) - 1)

m_truth <- moment_row(bt$M)
m_med   <- moment_row(bm$M)

fmt <- function(x) formatC(x, digits = 3, format = "g")
ratio <- function(a, b) if (b == 0) NA_real_ else a / b

moments_tbl <- tibble(
    Quantity   = c("E[M]", "Var[M]", "D = Var/E", "D - 1"),
    Truth      = c(fmt(m_truth$EM), fmt(m_truth$VarM), fmt(m_truth$D), fmt(m_truth$D_minus_1)),
    MEDICC2    = c(fmt(m_med$EM),   fmt(m_med$VarM),   fmt(m_med$D),   fmt(m_med$D_minus_1)),
    `M / T`    = c(fmt(ratio(m_med$EM,        m_truth$EM)),
                   fmt(ratio(m_med$VarM,      m_truth$VarM)),
                   "-",
                   fmt(ratio(m_med$D_minus_1, m_truth$D_minus_1)))
)

col_truth <- "#6BAED6"; col_med <- "#08306B"

p_branch <- ggplot(branch_df, aes(x = len, fill = source)) +
    geom_histogram(alpha = 0.5, position = "identity",
                   binwidth = max(1, floor(max(branch_df$len, na.rm = TRUE) / 30))) +
    geom_vline(xintercept = max(ec_truth$events),
               colour = col_truth, linetype = "dashed") +
    scale_fill_manual(values = c(truth = col_truth, MEDICC2 = col_med)) +
    scale_y_continuous(trans = "log1p",
                       breaks = c(0, 1, 3, 10, 30, 100)) +
    labs(title = "Branch length per edge",
         x = "events per edge", y = "count (log1p)",
         fill = NULL) +
    theme_bw(base_size = 11)

p_burden <- ggplot(burden_df, aes(x = M, fill = source)) +
    geom_histogram(alpha = 0.5, position = "identity",
                   binwidth = max(1, floor(max(burden_df$M, na.rm = TRUE) / 30))) +
    scale_fill_manual(values = c(truth = col_truth, MEDICC2 = col_med)) +
    labs(title = "Per-cell mutational burden M",
         x = "M (events on root->leaf path)", y = "cell count",
         fill = NULL) +
    theme_bw(base_size = 11)

p_moments <- gridExtra::tableGrob(
    moments_tbl,
    rows = NULL,
    theme = gridExtra::ttheme_minimal(base_size = 11)
)

title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(
        paste0("Parsimony-saturation diagnostics\n", sim_rel),
        fontface = "bold", size = 12
    )

fig <- cowplot::plot_grid(
    title_grob,
    p_branch,
    p_burden,
    p_moments,
    ncol = 1,
    rel_heights = c(0.08, 0.32, 0.32, 0.28)
)

out_path <- file.path(OUT_DIR, paste0("saturation__", sim_slug(sim_rel), ".png"))
ggsave(out_path, fig, width = 8, height = 12, dpi = 150)
message("wrote: ", out_path)
message("")
message("Moments summary:")
print(moments_tbl)
