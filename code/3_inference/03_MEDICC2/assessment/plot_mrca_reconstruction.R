# plot_mrca_reconstruction.R
#
# Answers three questions about MEDICC2's sample-MRCA CN reconstruction, using
# the per-sim stats cached by mrca_reconstruction_metrics.R:
#
#   (a) How often is the truth MRCA diploid? How often does MEDICC2 predict
#       diploid correctly?
#   (b) How often does MEDICC2 hallucinate CNAs at the MRCA, and how large are
#       those hallucinations?
#   (c) When the truth MRCA is non-diploid, how well does MEDICC2 reconstruct
#       it?
#
# All comparisons are per haplotype-bin (hbin): each MEDICC2 1 Mb bin
# contributes two positions (one per haplotype). "Diploid" means every hbin has
# CN = 1.
#
# Figures written to figures/3_inference/03_MEDICC2/MRCA_test/:
#   mrca_classification.png -- (a) sim-level 2x2 confusion matrix + per-scenario
#                              breakdown.
#   mrca_hallucinations.png -- (b) among truth-diploid MRCAs: how often does
#                              MEDICC2 hallucinate any alterations, and how
#                              many hbins/CN units are hallucinated?
#   mrca_altered_scatter.png -- (c) among truth-non-diploid MRCAs: bin-level
#                               truth vs MEDICC2 scatter + FP/FN counts.
#
# Usage:
#   Rscript code/3_inference/03_MEDICC2/assessment/plot_mrca_reconstruction.R

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(readr); library(ggplot2); library(cowplot)
})

args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE))
if (length(this_file) != 1L) stop("cannot resolve --file")
ROOT <- normalizePath(file.path(dirname(this_file), "..", "..", "..", ".."))

CSV_PATH <- file.path(ROOT, "data", "MEDICC2", "benchmark", "mrca_reconstruction.csv")
OUT_DIR  <- file.path(ROOT, "figures", "3_inference", "03_MEDICC2", "MRCA_test")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(CSV_PATH)) {
    stop("metrics CSV not found: ", CSV_PATH,
         "\nRun mrca_reconstruction_metrics.R first.")
}
df <- read_csv(CSV_PATH, show_col_types = FALSE)
message(sprintf("loaded %d rows from %s", nrow(df), CSV_PATH))

df <- df %>%
    mutate(
        scenario = factor(scenario,
                          levels = c("neutral", "selection_1", "selection_2")),
        truth_cls = ifelse(truth_is_diploid, "truth diploid", "truth altered"),
        med_cls   = ifelse(med_is_diploid,   "MEDICC2 diploid", "MEDICC2 altered"),
        truth_cls = factor(truth_cls, levels = c("truth diploid", "truth altered")),
        med_cls   = factor(med_cls,   levels = c("MEDICC2 diploid", "MEDICC2 altered"))
    )

scen_colours <- c(neutral     = "#66c2a5",
                  selection_1 = "#fc8d62",
                  selection_2 = "#8da0cb")

# ==========================================================================
# (a) Sim-level classification
# ==========================================================================
conf <- df %>%
    count(truth_cls, med_cls, name = "n") %>%
    complete(truth_cls, med_cls, fill = list(n = 0)) %>%
    mutate(
        cell_label = case_when(
            truth_cls == "truth diploid" & med_cls == "MEDICC2 diploid" ~ "TN\n(correct diploid)",
            truth_cls == "truth diploid" & med_cls == "MEDICC2 altered" ~ "FP\n(hallucination)",
            truth_cls == "truth altered" & med_cls == "MEDICC2 diploid" ~ "FN\n(missed)",
            truth_cls == "truth altered" & med_cls == "MEDICC2 altered" ~ "TP\n(both altered)"
        )
    )

p_a1 <- ggplot(conf, aes(x = med_cls, y = truth_cls, fill = n)) +
    geom_tile(colour = "grey30", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%s\nn = %d", cell_label, n)),
              size = 4, fontface = "bold", lineheight = 1.0) +
    scale_fill_gradient(low = "#f7fbff", high = "#2171b5", guide = "none") +
    scale_y_discrete(limits = rev) +
    labs(x = NULL, y = NULL,
         title = "(a) Sim-level MRCA classification: MEDICC2 vs truth",
         subtitle = sprintf("Each cell counts sims where the whole-MRCA state matches the given (truth, MEDICC2) pair. n = %d sims total. \"Altered\" = at least one hbin != diploid.", nrow(df))) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30", size = 9),
          axis.text = element_text(face = "bold"))

# Per-scenario stacked-bar breakdown
by_scen <- df %>%
    count(scenario, truth_cls, med_cls) %>%
    mutate(cls = case_when(
        truth_cls == "truth diploid" & med_cls == "MEDICC2 diploid" ~ "TN: correct diploid",
        truth_cls == "truth diploid" & med_cls == "MEDICC2 altered" ~ "FP: hallucination",
        truth_cls == "truth altered" & med_cls == "MEDICC2 diploid" ~ "FN: missed",
        truth_cls == "truth altered" & med_cls == "MEDICC2 altered" ~ "TP: both altered"
    ),
    cls = factor(cls, levels = c("TN: correct diploid", "TP: both altered",
                                 "FP: hallucination", "FN: missed")))

cls_cols <- c("TN: correct diploid" = "#a6cee3",
              "TP: both altered"    = "#33a02c",
              "FP: hallucination"   = "#e31a1c",
              "FN: missed"          = "#ff7f00")

p_a2 <- ggplot(by_scen, aes(x = scenario, y = n, fill = cls)) +
    geom_col(position = position_stack(), colour = "grey20", linewidth = 0.2) +
    geom_text(aes(label = ifelse(n > 0, n, "")),
              position = position_stack(vjust = 0.5), size = 3.2, colour = "white",
              fontface = "bold") +
    scale_fill_manual(values = cls_cols, name = NULL, drop = FALSE) +
    labs(x = NULL, y = "Number of sims",
         title = "Per-scenario breakdown") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "right",
          legend.key.height = unit(0.8, "lines"))

fig_a <- plot_grid(p_a1, p_a2, ncol = 2, rel_widths = c(1, 1.05), align = "h", axis = "tb")
out_a <- file.path(OUT_DIR, "mrca_classification.png")
ggsave(out_a, fig_a, width = 13, height = 5.2, dpi = 130)
message("wrote: ", out_a)

# ==========================================================================
# (b) Hallucinations among truth-diploid MRCAs
# ==========================================================================
truth_dip <- df %>% filter(truth_is_diploid)
n_td <- nrow(truth_dip)
n_hal <- sum(truth_dip$halluc_hbins > 0)

# Left panel: fraction with any hallucination
frac_df <- truth_dip %>%
    mutate(any_halluc = halluc_hbins > 0) %>%
    count(scenario, any_halluc) %>%
    complete(scenario, any_halluc = c(FALSE, TRUE), fill = list(n = 0)) %>%
    mutate(halluc_lbl = ifelse(any_halluc, "hallucinated", "clean diploid"),
           halluc_lbl = factor(halluc_lbl, levels = c("clean diploid", "hallucinated")))

p_b1 <- ggplot(frac_df, aes(x = scenario, y = n, fill = halluc_lbl)) +
    geom_col(colour = "grey20", linewidth = 0.2) +
    geom_text(aes(label = ifelse(n > 0, n, "")),
              position = position_stack(vjust = 0.5),
              size = 3.2, colour = "white", fontface = "bold") +
    scale_fill_manual(values = c("clean diploid" = "#a6cee3",
                                 "hallucinated"  = "#e31a1c"),
                      name = NULL) +
    labs(x = NULL, y = "Number of sims (truth = diploid)",
         title = sprintf("(b1) Among truth-diploid MRCAs (n=%d), does MEDICC2 hallucinate?", n_td),
         subtitle = sprintf("%d / %d truth-diploid sims have at least one MEDICC2-altered hbin.", n_hal, n_td)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30", size = 9),
          legend.position = "right")

# Right panel: severity of the hallucinations that DO occur
hal_only <- truth_dip %>% filter(halluc_hbins > 0)
p_b2 <- if (nrow(hal_only) == 0L) {
    ggplot() + annotate("text", x = 0, y = 0, size = 6,
                        label = "No hallucinations occurred.\nAll truth-diploid MRCAs\nreconstructed cleanly.") +
        theme_void() +
        labs(title = "(b2) Hallucination severity") +
        theme(plot.title = element_text(face = "bold"))
} else {
    ggplot(hal_only, aes(x = halluc_hbins, y = halluc_units, colour = scenario)) +
        geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dashed") +
        geom_point(size = 3.5, alpha = 0.9) +
        geom_text(aes(label = sprintf("%s / %s\nN=%s d=%s", timing, sim_id, N,
                                      ifelse(is.na(d), "-", as.character(d)))),
                  size = 2.8, colour = "grey20",
                  hjust = -0.1, vjust = -0.4, lineheight = 0.9) +
        scale_colour_manual(values = scen_colours, name = NULL, drop = FALSE) +
        expand_limits(x = 0, y = 0) +
        scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
        labs(x = "Hallucinated hbins  (# positions where MEDICC2 says altered, truth says diploid)",
             y = "Hallucinated CN units  (Sum |MEDICC2 - 1| over those hbins)",
             title = "(b2) Hallucination severity for the sims that do hallucinate",
             subtitle = "y = x (dashed) => hallucinations are all magnitude-1 (each altered hbin off by exactly 1 CN unit).") +
        theme_bw(base_size = 11) +
        theme(panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(colour = "grey30", size = 9))
}

fig_b <- plot_grid(p_b1, p_b2, ncol = 2, rel_widths = c(1, 1.3), align = "h", axis = "tb")
out_b <- file.path(OUT_DIR, "mrca_hallucinations.png")
ggsave(out_b, fig_b, width = 14, height = 5.2, dpi = 130)
message("wrote: ", out_b)

# ==========================================================================
# (c) Non-diploid MRCA reconstruction quality
# ==========================================================================
alt <- df %>% filter(!truth_is_diploid)
n_alt <- nrow(alt)

# Scatter: truth vs MEDICC2 altered hbins
p_c1 <- ggplot(alt, aes(x = truth_altered_hbins, y = med_altered_hbins,
                        colour = scenario)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dashed") +
    geom_point(size = 2.4, alpha = 0.8) +
    scale_colour_manual(values = scen_colours, name = NULL, drop = FALSE) +
    expand_limits(x = 0, y = 0) +
    coord_fixed() +
    labs(x = "Truth altered hbins",
         y = "MEDICC2 altered hbins",
         title = sprintf("(c1) Non-diploid MRCAs (n=%d): total altered-hbin count", n_alt),
         subtitle = "Points on y = x (dashed) mean MEDICC2 flagged the same number of hbins as truth.") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30", size = 9))

# FP / FN barcode: how many of MEDICC2's altered hbins were spurious, how many
# truth-altered hbins did MEDICC2 miss?
err_long <- alt %>%
    select(sim_rel, scenario, fp_hbins, fn_hbins, miscall_hbins) %>%
    pivot_longer(c(fp_hbins, fn_hbins, miscall_hbins),
                 names_to = "kind", values_to = "hbins") %>%
    mutate(kind = recode(kind,
                         fp_hbins      = "FP (hallucinated)",
                         fn_hbins      = "FN (missed)",
                         miscall_hbins = "Miscall (both altered, wrong CN)"),
           kind = factor(kind,
                         levels = c("FP (hallucinated)",
                                    "FN (missed)",
                                    "Miscall (both altered, wrong CN)")))

n_perfect_alt <- alt %>% filter(fp_hbins == 0, fn_hbins == 0, miscall_hbins == 0) %>% nrow()

p_c2 <- ggplot(err_long, aes(x = kind, y = hbins, colour = scenario)) +
    geom_hline(yintercept = 0, colour = "grey60", linetype = "dashed") +
    geom_jitter(width = 0.18, height = 0, size = 2.2, alpha = 0.75) +
    scale_colour_manual(values = scen_colours, name = NULL, drop = FALSE) +
    labs(x = NULL, y = "Error hbins per sim",
         title = "(c2) Where the reconstruction disagrees",
         subtitle = sprintf("%d / %d non-diploid sims are perfect bin-by-bin (all three error kinds = 0).",
                            n_perfect_alt, n_alt)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30", size = 9),
          axis.text.x = element_text(size = 9))

fig_c <- plot_grid(p_c1, p_c2, ncol = 2, rel_widths = c(1, 1.1), align = "h", axis = "tb")
out_c <- file.path(OUT_DIR, "mrca_altered_scatter.png")
ggsave(out_c, fig_c, width = 13, height = 5.2, dpi = 130)
message("wrote: ", out_c)

# ==========================================================================
# Console summary
# ==========================================================================
cat("\n---- Summary ----\n")
cat(sprintf("Total sims: %d\n", nrow(df)))
cat(sprintf("Truth MRCA diploid:   %d (%.1f%%)\n",
            sum(df$truth_is_diploid), 100 * mean(df$truth_is_diploid)))
cat(sprintf("MEDICC2 MRCA diploid: %d (%.1f%%)\n",
            sum(df$med_is_diploid), 100 * mean(df$med_is_diploid)))
cat(sprintf("Exact bin-by-bin match: %d (%.1f%%)\n",
            sum(df$exact_match), 100 * mean(df$exact_match)))
cat(sprintf("\n(a) Sim-level confusion:\n"))
print(with(df, table(truth = truth_cls, medicc = med_cls)))
cat(sprintf("\n(b) Among %d truth-diploid sims, %d had any hallucination.\n", n_td, n_hal))
if (nrow(hal_only) > 0L) {
    cat("Hallucinating sims:\n")
    print(hal_only %>% select(sim_rel, halluc_hbins, halluc_units,
                              halluc_gain_hbins, halluc_loss_hbins))
}
cat(sprintf("\n(c) Among %d truth-non-diploid sims, %d are perfect bin-by-bin.\n",
            n_alt, n_perfect_alt))
imperfect_alt <- alt %>% filter(fp_hbins > 0 | fn_hbins > 0 | miscall_hbins > 0)
if (nrow(imperfect_alt) > 0L) {
    cat("Imperfect non-diploid reconstructions:\n")
    print(imperfect_alt %>%
              select(sim_rel, truth_altered_hbins, med_altered_hbins,
                     tp_hbins, fp_hbins, fn_hbins, miscall_hbins))
}
