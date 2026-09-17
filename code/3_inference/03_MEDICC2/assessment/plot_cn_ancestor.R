# plot_cn_ancestor.R
#
# Aggregate figures for MEDICC2 internal-node CN-profile reconstruction. Reads
# both tiers of output from cn_ancestor_metrics.R:
#   - data/MEDICC2/benchmark/cn_ancestor_summary.csv  (one row per sim)
#   - data/MEDICC2/benchmark/cn_ancestor/<slug>.csv   (per-sim per-node metrics)
#
# Also joins truth_EM from data/MEDICC2/benchmark/aggregate_metrics.csv so the
# by-scenario scatter can share the x-axis of the RF aggregate figures.
#
# Figures written to figures/3_inference/03_MEDICC2/MRCA_test/:
#   (1) cn_metric_comparison.png       -- boxplots of all 6 sample-MRCA metrics
#                                         per scenario, so the three
#                                         normalisation choices can be compared
#                                         at a glance.
#   (2) cn_sample_mrca_neutral.png     -- neutral-only bar grid (cols = d,
#                                         bars per timing x N) using the
#                                         baseline-relative Hamming as the
#                                         headline metric. Mirrors
#                                         neutral_rf_grid.png.
#   (3) cn_sample_mrca_by_scenario.png -- 3-panel scatter: baseline-relative
#                                         Hamming vs truth E[M], per scenario.
#                                         Mirrors rf_by_scenario.png.
#   (4) cn_stratification.png          -- clade-size stratification, 2 rows
#                                         (strategy A / strategy B) x 3 cols
#                                         (scenario). Direct test of "worse
#                                         closer to root".
#
# Usage:
#   Rscript code/3_inference/03_MEDICC2/assessment/plot_cn_ancestor.R

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(readr); library(ggplot2); library(cowplot)
})

args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE))
if (length(this_file) != 1L) stop("cannot resolve --file")
ROOT <- normalizePath(file.path(dirname(this_file), "..", "..", "..", ".."))

SUMMARY_CSV   <- file.path(ROOT, "data", "MEDICC2", "benchmark", "cn_ancestor_summary.csv")
PER_SIM_DIR   <- file.path(ROOT, "data", "MEDICC2", "benchmark", "cn_ancestor")
AGGREGATE_CSV <- file.path(ROOT, "data", "MEDICC2", "benchmark", "aggregate_metrics.csv")
OUT_DIR <- file.path(ROOT, "figures", "3_inference", "03_MEDICC2",
                     "MRCA_test")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(SUMMARY_CSV)) {
    stop("summary CSV not found: ", SUMMARY_CSV,
         "\nRun cn_ancestor_metrics.R first.")
}
summ <- read_csv(SUMMARY_CSV, show_col_types = FALSE)
message(sprintf("loaded %d sims from %s", nrow(summ), SUMMARY_CSV))

agg <- if (file.exists(AGGREGATE_CSV)) {
    read_csv(AGGREGATE_CSV, show_col_types = FALSE) %>%
        select(sim_rel, truth_EM)
} else {
    warning("aggregate_metrics.csv missing -- truth_EM will be NA")
    tibble(sim_rel = character(), truth_EM = numeric())
}
summ <- summ %>% left_join(agg, by = "sim_rel")

# ---- Shared palettes (mirror plot_aggregate_metrics.R) ---------------------
timing_colours <- c(
    deterministic = "#4682B4",
    gamma         = "#9370DB",
    markov        = "#008080"
)
shape_map <- c("0.0" = 17, "0.5" = 16, "0.9" = 15, "n/a" = 8)

summ <- summ %>%
    mutate(
        d_lab = case_when(
            is.na(d) ~ "n/a",
            TRUE     ~ sprintf("%.1f", d)
        ),
        d_lab = factor(d_lab, levels = c("0.0", "0.5", "0.9", "n/a")),
        N_group = factor(ifelse(N < 5000, "N ~ 1k", "N ~ 10k+"),
                         levels = c("N ~ 1k", "N ~ 10k+")),
        timing = factor(timing, levels = c("deterministic", "gamma", "markov"))
    )

# ==========================================================================
# Figure 1: sample-MRCA metric comparison (all 6 metrics side-by-side)
# ==========================================================================
metric_long <- summ %>%
    select(sim_rel, scenario, timing,
           sample_mrca_hamming_raw, sample_mrca_hamming_signal, sample_mrca_hamming_baseline,
           sample_mrca_sumabs_raw,  sample_mrca_sumabs_signal,  sample_mrca_sumabs_baseline) %>%
    pivot_longer(starts_with("sample_mrca_"),
                 names_to = "metric", values_to = "value") %>%
    mutate(
        base = ifelse(grepl("hamming", metric), "Hamming", "Sum-|Delta|"),
        norm = case_when(
            grepl("_raw$",      metric) ~ "raw",
            grepl("_signal$",   metric) ~ "signal-restricted",
            grepl("_baseline$", metric) ~ "baseline-relative",
            TRUE ~ NA_character_
        ),
        base = factor(base, levels = c("Hamming", "Sum-|Delta|")),
        norm = factor(norm, levels = c("raw", "signal-restricted", "baseline-relative"))
    )

p1 <- ggplot(metric_long, aes(x = scenario, y = value, fill = scenario)) +
    geom_boxplot(outlier.size = 0.6, alpha = 0.7, linewidth = 0.3) +
    geom_jitter(width = 0.15, size = 0.4, alpha = 0.5, colour = "grey20") +
    facet_grid(base ~ norm, scales = "free_y", switch = "y") +
    scale_fill_manual(values = c(neutral = "#66c2a5",
                                 selection_1 = "#fc8d62",
                                 selection_2 = "#8da0cb"),
                      name = NULL) +
    labs(x = NULL, y = NULL,
         title = "Sample-MRCA CN reconstruction -- comparison of six metrics",
         subtitle = "Rows: base distance. Cols: normalisation. Each dot = one sim. baseline-relative (right) is the recommended headline metric: 1 = perfect, 0 = no better than 'all diploid', negative = worse.") +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"),
          strip.placement = "outside",
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30"),
          axis.text.x = element_text(angle = 20, hjust = 1))

out1 <- file.path(OUT_DIR, "cn_metric_comparison.png")
ggsave(out1, p1, width = 12, height = 7, dpi = 130)
message("wrote: ", out1)

# ==========================================================================
# Figure 2: Neutral-only sample-MRCA bar grid (baseline-relative Hamming)
# ==========================================================================
neutral <- summ %>% filter(scenario == "neutral")

det <- neutral %>% filter(timing == "deterministic")
non_det <- neutral %>% filter(timing != "deterministic")
det_expanded <- bind_rows(det %>% mutate(d = 0.0),
                          det %>% mutate(d = 0.5),
                          det %>% mutate(d = 0.9))
grid_df <- bind_rows(non_det, det_expanded) %>%
    mutate(d_lab = factor(sprintf("d = %.1f", d),
                          levels = c("d = 0.0", "d = 0.5", "d = 0.9")))

p2 <- ggplot(grid_df,
             aes(x = timing, y = sample_mrca_hamming_baseline,
                 fill = timing, alpha = N_group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             colour = "grey20", linewidth = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    facet_wrap(~ d_lab, nrow = 1) +
    scale_fill_manual(values = timing_colours, name = "Timing") +
    scale_alpha_manual(values = c("N ~ 1k" = 0.55, "N ~ 10k+" = 1.0),
                       name = "Population size") +
    coord_cartesian(ylim = c(min(0, min(grid_df$sample_mrca_hamming_baseline, na.rm = TRUE)),
                             1.05)) +
    labs(x = NULL, y = "Sample-MRCA baseline-relative Hamming",
         title = "Neutral simulations - MEDICC2 sample-MRCA CN reconstruction",
         subtitle = "1 = perfect reconstruction (dashed), 0 = no better than 'all diploid'. Bars = timing x N. Deterministic bars repeat across d columns since deterministic sims have no d parameter.") +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30"))

out2 <- file.path(OUT_DIR, "cn_sample_mrca_neutral.png")
ggsave(out2, p2, width = 11, height = 4.5, dpi = 130)
message("wrote: ", out2)

# ==========================================================================
# Figure 3: Sample-MRCA by scenario -- baseline-relative Hamming vs truth_EM
# ==========================================================================
scen_panel <- function(sub, color_kind = c("timing", "s"), title) {
    color_kind <- match.arg(color_kind)
    p <- ggplot(sub, aes(x = truth_EM, y = sample_mrca_hamming_baseline)) +
        geom_hline(yintercept = 1, colour = "grey60", linetype = "dashed") +
        geom_hline(yintercept = 0, colour = "grey80", linetype = "dotted")
    if (color_kind == "s" && "s" %in% names(sub) && any(!is.na(sub$s))) {
        p <- p + geom_line(
            data = sub %>% arrange(s),
            aes(group = interaction(timing, N, d_lab), colour = s),
            alpha = 0.4, linewidth = 0.4
        )
    }
    p <- p + geom_point(
        aes(colour = if (color_kind == "timing") timing else s,
            shape = d_lab, size = N),
        alpha = 0.85
    )
    if (color_kind == "timing") {
        p <- p + scale_colour_manual(values = timing_colours, name = "Timing", na.value = "grey40")
    } else {
        p <- p + scale_colour_viridis_c(name = "Selection s")
    }
    p +
        scale_shape_manual(values = shape_map, name = "Death rate d", drop = FALSE) +
        scale_size_continuous(name = "N", breaks = c(1000, 10000, 16384),
                              range = c(1.6, 3.6), trans = "log10") +
        coord_cartesian(ylim = c(min(-0.05, min(summ$sample_mrca_hamming_baseline, na.rm = TRUE)),
                                 1.05)) +
        labs(x = "Truth E[M] (mean per-cell root->leaf events)",
             y = "Baseline-relative Hamming (sample-MRCA)",
             title = title) +
        theme_bw(base_size = 10) +
        theme(panel.grid.minor = element_blank(),
              legend.position = "right",
              legend.key.height = unit(0.6, "lines"),
              plot.title = element_text(face = "bold"),
              aspect.ratio = 1) +
        guides(shape = "none", size = "none")
}

scen_order <- c("neutral", "selection_1", "selection_2")
panels <- lapply(scen_order, function(sc) {
    sub <- summ %>% filter(scenario == sc)
    ck <- if (sc == "neutral") "timing" else "s"
    title <- switch(sc,
                    neutral     = sprintf("Neutral (n=%d)",     nrow(sub)),
                    selection_1 = sprintf("Selection 1 (n=%d)", nrow(sub)),
                    selection_2 = sprintf("Selection 2 (n=%d)", nrow(sub)))
    scen_panel(sub, ck, title)
})

title_grob <- cowplot::ggdraw() + cowplot::draw_label(
    "Sample-MRCA CN reconstruction vs tree depth",
    fontface = "bold", size = 13, hjust = 0.5
)
subtitle_grob <- cowplot::ggdraw() + cowplot::draw_label(
    "1 = perfect (dashed); 0 = no better than 'all diploid' (dotted). Deeper truth trees (large E[M]) should be harder to reconstruct at the root.",
    size = 10, hjust = 0.5, colour = "grey30"
)

legend_plot <- scen_panel(summ, "timing", "") +
    guides(colour = "none") +
    theme(legend.position = "bottom", legend.box = "horizontal")
shared_legend <- cowplot::get_legend(legend_plot)

body <- cowplot::plot_grid(plotlist = panels, ncol = 3, align = "h")
fig3 <- cowplot::plot_grid(title_grob, subtitle_grob, body, shared_legend,
                           ncol = 1, rel_heights = c(0.05, 0.04, 1, 0.08))
out3 <- file.path(OUT_DIR, "cn_sample_mrca_by_scenario.png")
ggsave(out3, fig3, width = 15, height = 5.5, dpi = 130)
message("wrote: ", out3)

# ==========================================================================
# Figure 4: Stratification by clade size (pooled per-node CSVs)
# ==========================================================================
per_sim_files <- list.files(PER_SIM_DIR, pattern = "\\.csv$", full.names = TRUE)
if (length(per_sim_files) == 0L) {
    warning("no per-sim CSVs under ", PER_SIM_DIR, " -- skipping stratification figure")
} else {
    # Slug -> sim_rel via the summary table
    slug_of <- function(sim_rel) {
        parts <- strsplit(sim_rel, "/", fixed = TRUE)[[1]]
        paste(tail(parts, 2), collapse = "__")
    }
    slug_to_rel <- setNames(summ$sim_rel, vapply(summ$sim_rel, slug_of, character(1)))

    read_sim_file <- function(path) {
        slug <- sub("\\.csv$", "", basename(path))
        sim_rel <- slug_to_rel[[slug]]
        if (is.null(sim_rel)) return(NULL)
        d <- read_csv(path, show_col_types = FALSE,
                      col_types = cols(.default = "d", strategy = "c",
                                       leaf_i = "c", leaf_j = "c",
                                       truth_mrca_name = "c", medicc_mrca_name = "c"))
        d$sim_rel <- sim_rel
        d
    }

    all_rows <- bind_rows(lapply(per_sim_files, read_sim_file)) %>%
        left_join(summ %>% select(sim_rel, scenario, timing, N, d, s),
                  by = "sim_rel") %>%
        mutate(scenario = factor(scenario, levels = scen_order),
               timing   = factor(timing, levels = c("deterministic", "gamma", "markov")))

    message(sprintf("pooled %d per-sim CSVs (%d rows)",
                    length(per_sim_files), nrow(all_rows)))

    # For strategy B, deduplicate per sim to unique truth_mrca so nodes
    # subtending many pairs aren't over-plotted. Strategy A rows are already
    # unique per matched clade.
    plot_df <- all_rows %>%
        filter(strategy %in% c("A", "B")) %>%
        mutate(node_key = paste(strategy, sim_rel, truth_mrca_name, sep = "|")) %>%
        distinct(node_key, .keep_all = TRUE) %>%
        select(-node_key)

    strat_labels <- c(A = "Strategy A (matched clades)",
                      B = "Strategy B (pair-MRCA, deduplicated)")
    plot_df$strategy_lab <- factor(strat_labels[plot_df$strategy],
                                   levels = unname(strat_labels))

    p4 <- ggplot(plot_df,
                 aes(x = clade_size_truth, y = hamming_baseline,
                     colour = timing)) +
        geom_hline(yintercept = 1, colour = "grey60", linetype = "dashed") +
        geom_hline(yintercept = 0, colour = "grey80", linetype = "dotted") +
        geom_point(alpha = 0.1, size = 0.5) +
        geom_smooth(aes(group = timing), method = "loess", se = TRUE,
                    linewidth = 0.7, alpha = 0.15) +
        facet_grid(strategy_lab ~ scenario) +
        scale_x_log10(breaks = c(2, 5, 10, 30, 100, 500)) +
        scale_colour_manual(values = timing_colours, name = "Timing",
                            na.value = "grey40") +
        coord_cartesian(ylim = c(-0.2, 1.05)) +
        labs(x = "Truth-MRCA clade size (log scale)",
             y = "Baseline-relative Hamming",
             title = "CN reconstruction quality vs tree depth (clade-size stratification)",
             subtitle = "Rows: matching strategy. Cols: scenario. Larger clade size = closer to root. 1 = perfect (dashed); 0 = 'all diploid' baseline (dotted).") +
        theme_bw(base_size = 10) +
        theme(strip.background = element_rect(fill = "grey95", colour = NA),
              strip.text = element_text(face = "bold"),
              panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold"),
              plot.subtitle = element_text(colour = "grey30"),
              legend.position = "bottom")

    out4 <- file.path(OUT_DIR, "cn_stratification.png")
    ggsave(out4, p4, width = 12, height = 7, dpi = 130)
    message("wrote: ", out4)
}
