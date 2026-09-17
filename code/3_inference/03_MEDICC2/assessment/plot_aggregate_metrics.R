# plot_aggregate_metrics.R
#
# Aggregate figures over all MEDICC2 sims cached by aggregate_metrics.R.
#
# Parsimony-saturation (-> figures/.../parsimony_saturation/aggregate/):
#   (1) branch_length_agreement.png : truth vs MEDICC2 mean branch length per edge
#   (2) burden_agreement.png        : truth vs MEDICC2 mean per-cell burden E[M]
#   (3) cophenetic_slope.png        : slope of MEDICC2 cophenetic distances on
#                                     truth, plotted against truth E[M]
#   (4) neutral_summary_grid.png    : neutral-only 3x3 grid (rows = stat,
#                                     cols = d), bars per timing x N.
#
# Phylogeny topology (-> figures/.../phylogenies/aggregate/):
#   (5) neutral_rf_grid.png         : neutral-only normalised RF, cols = d,
#                                     bars per timing x N.
#   (6) rf_by_scenario.png          : 3 panels (neutral / sel1 / sel2), y =
#                                     normalised RF vs truth E[M].
#
# Each 3-panel figure is faceted by scenario (neutral | selection_1 | selection_2)
# so each panel can use a colour scale that carries information for that scenario:
#   - neutral       : colour = timing model (discrete)
#   - selection_1/2 : colour = selection strength s (viridis)
# Point shape encodes death rate d; point size encodes population size N.
# Sel1/sel2 points within a (timing, N, d) group are connected by a line
# ordered by s, so the saturation trajectory as selection strengthens is
# visible.
#
# Usage:
#   Rscript code/3_inference/03_MEDICC2/assessment/plot_aggregate_metrics.R

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(readr); library(ggplot2); library(cowplot)
})

args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE))
if (length(this_file) != 1L) stop("cannot resolve --file")
ROOT <- normalizePath(file.path(dirname(this_file), "..", "..", "..", ".."))

source(file.path(ROOT, "code", "X_helpers", "plotting_functions.R"))

CSV_PATH <- file.path(ROOT, "data", "MEDICC2", "benchmark", "aggregate_metrics.csv")
OUT_DIR  <- file.path(ROOT, "figures", "3_inference", "03_MEDICC2",
                      "parsimony_saturation", "aggregate")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(CSV_PATH)) {
    stop("aggregate CSV not found: ", CSV_PATH,
         "\nRun code/3_inference/03_MEDICC2/assessment/aggregate_metrics.R first.")
}
df <- read_csv(CSV_PATH, show_col_types = FALSE)
message(sprintf("loaded %d rows from %s", nrow(df), CSV_PATH))

# Death rate as an ordered discrete for shape mapping. Deterministic sims have
# NA for d -- give them their own bucket so they still plot.
df <- df %>%
    mutate(
        d_lab = case_when(
            is.na(d) ~ "n/a",
            TRUE     ~ sprintf("%.1f", d)
        ),
        d_lab = factor(d_lab, levels = c("0.0", "0.5", "0.9", "n/a"))
    )

# Palette mirroring the Julia code in plotting_functions.jl / .R.
timing_colours <- c(
    deterministic = "#4682B4",
    gamma         = "#9370DB",
    markov        = "#008080"
)
shape_map <- c("0.0" = 17, "0.5" = 16, "0.9" = 15, "n/a" = 8)

# ---------------------------------------------------------------------------
# Panel builder. `color_kind` controls whether the panel uses discrete timing
# colours (neutral) or a viridis-continuous s scale (selection).
# ---------------------------------------------------------------------------
panel <- function(sub, xvar, yvar, xlab, ylab, title,
                  color_kind = c("timing", "s"),
                  diagonal = TRUE, ref_y = 1) {
    color_kind <- match.arg(color_kind)
    p <- ggplot(sub, aes(x = .data[[xvar]], y = .data[[yvar]]))

    if (diagonal) {
        p <- p + geom_abline(slope = 1, intercept = 0,
                             colour = "grey60", linetype = "dashed")
    } else {
        # Non-diagonal reference at y = ref_y (perfect = 1 for cophenetic
        # slope, 0 for normalised RF; caller decides).
        p <- p + geom_hline(yintercept = ref_y, colour = "grey60", linetype = "dashed")
    }

    # For selection panels, draw a line per (timing, N, d) group ordered by s
    # so the saturation trajectory is visible. Neutrals have no s -> no line.
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
        p <- p + scale_colour_manual(values = timing_colours, name = "Timing",
                                     na.value = "grey40")
    } else {
        p <- p + scale_colour_viridis_c(name = "Selection s", option = "viridis")
    }

    p +
        scale_shape_manual(values = shape_map, name = "Death rate d",
                           drop = FALSE) +
        scale_size_continuous(name = "N",
                              breaks = c(1000, 10000, 16384),
                              range = c(1.6, 3.6),
                              trans = "log10") +
        labs(x = xlab, y = ylab, title = title) +
        theme_bw(base_size = 10) +
        theme(panel.grid.minor = element_blank(),
              legend.position  = "right",
              legend.key.height = unit(0.6, "lines"),
              legend.box.spacing = unit(0.2, "lines"),
              plot.title = element_text(face = "bold"))
}

# ---------------------------------------------------------------------------
# Assemble one figure (3 panels: neutral, sel1, sel2). Each panel scales to
# its own data so low-burden selection panels are readable even when a deep
# neutral sim (e.g. N=16384) pushes the neutral axis to hundreds. Within
# each panel coord_fixed keeps y=x an honest 45 deg reference.
# ---------------------------------------------------------------------------
build_figure <- function(xvar, yvar, xlab, ylab, big_title,
                         subtitle, diagonal = TRUE,
                         ref_y = 1, ylim_override = NULL) {
    scen_order <- c("neutral", "selection_1", "selection_2")

    # For non-diagonal figures use one common y-range so vertical position is
    # directly comparable across panels. `ylim_override` wins when set (e.g.
    # RF is bounded to [0, 1] by definition); otherwise derive from data.
    common_yrange <- if (!diagonal) {
        if (!is.null(ylim_override)) ylim_override
        else {
            yr <- range(df[[yvar]], na.rm = TRUE, finite = TRUE)
            yr + c(-1, 1) * max(diff(yr) * 0.05, 0.01)
        }
    } else NULL

    panels <- lapply(scen_order, function(sc) {
        sub <- df %>% filter(scenario == sc)
        color_kind <- if (sc == "neutral") "timing" else "s"
        title <- switch(sc,
                        neutral       = sprintf("Neutral (n=%d)",       nrow(sub)),
                        selection_1   = sprintf("Selection 1 (n=%d)",   nrow(sub)),
                        selection_2   = sprintf("Selection 2 (n=%d)",   nrow(sub)))
        p <- panel(sub, xvar, yvar, xlab, ylab, title,
                   color_kind = color_kind, diagonal = diagonal, ref_y = ref_y)
        if (diagonal) {
            xy <- c(sub[[xvar]], sub[[yvar]])
            rng <- range(xy, na.rm = TRUE, finite = TRUE)
            pad <- max(diff(rng) * 0.05, 0.5)
            lims <- c(max(0, rng[1] - pad), rng[2] + pad)
            p <- p + coord_fixed(xlim = lims, ylim = lims)
        } else {
            xr <- range(sub[[xvar]], na.rm = TRUE, finite = TRUE)
            p <- p + coord_cartesian(
                xlim = xr + c(-1, 1) * max(diff(xr) * 0.05, 0.5),
                ylim = common_yrange
            ) + theme(aspect.ratio = 1)  # keep the three panels square
        }
        # Hide the shape+size legends on all panels; they'll be shown once at
        # the bottom of the figure. The color legend stays per-panel because
        # its meaning differs (timing for neutral, s for selection).
        p + guides(shape = "none", size = "none")
    })

    title_grob <- cowplot::ggdraw() + cowplot::draw_label(
        big_title, fontface = "bold", size = 13, hjust = 0.5
    )
    subtitle_grob <- cowplot::ggdraw() + cowplot::draw_label(
        subtitle, size = 10, hjust = 0.5, colour = "grey30"
    )

    # Extract a shape+size legend from a dummy plot that shows both.
    legend_plot <- panel(df, xvar, yvar, xlab, ylab, "",
                         color_kind = "timing", diagonal = diagonal,
                         ref_y = ref_y) +
        guides(colour = "none") +
        theme(legend.position = "bottom", legend.box = "horizontal")
    shared_legend <- cowplot::get_legend(legend_plot)

    body <- cowplot::plot_grid(plotlist = panels, ncol = 3, align = "h")
    cowplot::plot_grid(title_grob, subtitle_grob, body, shared_legend,
                       ncol = 1, rel_heights = c(0.05, 0.04, 1, 0.08))
}

# ---------------------------------------------------------------------------
# Emit the three figures
# ---------------------------------------------------------------------------
fig1 <- build_figure(
    xvar = "truth_bl_mean", yvar = "med_bl_mean",
    xlab = "Truth: mean events per collapsed edge",
    ylab = "MEDICC2: mean events per edge",
    big_title = "Mean branch length per edge -- truth vs MEDICC2",
    subtitle  = "Points below the y=x diagonal indicate parsimony compressing multiple truth events into a single MEDICC2 edge.",
    diagonal  = TRUE
)
out1 <- file.path(OUT_DIR, "branch_length_agreement.png")
ggsave(out1, fig1, width = 15, height = 5.5, dpi = 130)
message("wrote: ", out1)

fig2 <- build_figure(
    xvar = "truth_EM", yvar = "med_EM",
    xlab = "Truth E[M] (mean per-cell root->leaf events)",
    ylab = "MEDICC2 E[M]",
    big_title = "Mean per-cell mutational burden E[M] -- truth vs MEDICC2",
    subtitle  = "Points below y=x mean MEDICC2 underestimates burden along the root->leaf path (parsimony-saturation).",
    diagonal  = TRUE
)
out2 <- file.path(OUT_DIR, "burden_agreement.png")
ggsave(out2, fig2, width = 15, height = 5.5, dpi = 130)
message("wrote: ", out2)

fig3 <- build_figure(
    xvar = "truth_EM", yvar = "coph_slope",
    xlab = "Truth E[M] (mean per-cell root->leaf events)",
    ylab = "Cophenetic slope (MEDICC2 ~ truth)",
    big_title = "Cophenetic pairwise-distance slope (MEDICC2 on truth)",
    subtitle  = "Slope < 1 means MEDICC2's pairwise cell-cell distances are compressed relative to truth.",
    diagonal  = FALSE
)
out3 <- file.path(OUT_DIR, "cophenetic_slope.png")
ggsave(out3, fig3, width = 15, height = 5.5, dpi = 130)
message("wrote: ", out3)

# ---------------------------------------------------------------------------
# Neutral-only summary grid: rows = statistic, cols = death rate d.
# Within each cell, one bar per timing model x N-size combo. Ratios are
# MEDICC2 / truth so the "perfect" reference is y = 1 (dashed).
# Deterministic sims have no d parameter, so their bars are repeated in each
# d column to keep the timing comparison meaningful across the row.
# ---------------------------------------------------------------------------
neutral <- df %>% filter(scenario == "neutral") %>%
    mutate(
        em_ratio  = med_EM / truth_EM,
        bl_ratio  = med_bl_mean / truth_bl_mean,
        N_group   = factor(ifelse(N < 5000, "N ~ 1k", "N ~ 10k+"),
                           levels = c("N ~ 1k", "N ~ 10k+"))
    )

neutral_long <- neutral %>%
    select(timing, N, N_group, d, coph_slope, em_ratio, bl_ratio) %>%
    pivot_longer(c(coph_slope, em_ratio, bl_ratio),
                 names_to = "stat", values_to = "value") %>%
    mutate(
        stat = factor(stat,
                      levels = c("coph_slope", "em_ratio", "bl_ratio"),
                      labels = c("Cophenetic slope",
                                 "MEDICC2 / truth : E[M]",
                                 "MEDICC2 / truth : mean events per edge")),
        timing = factor(timing, levels = c("deterministic", "gamma", "markov"))
    )

det <- neutral_long %>% filter(timing == "deterministic")
non_det <- neutral_long %>% filter(timing != "deterministic")
det_expanded <- bind_rows(
    det %>% mutate(d = 0.0),
    det %>% mutate(d = 0.5),
    det %>% mutate(d = 0.9)
)
grid_df <- bind_rows(non_det, det_expanded) %>%
    mutate(d_lab = factor(sprintf("d = %.1f", d),
                          levels = c("d = 0.0", "d = 0.5", "d = 0.9")))

p_neutral_grid <- ggplot(grid_df,
                         aes(x = timing, y = value,
                             fill = timing, alpha = N_group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             colour = "grey20", linewidth = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    facet_grid(stat ~ d_lab, scales = "free_y", switch = "y") +
    scale_fill_manual(values = timing_colours, name = "Timing") +
    scale_alpha_manual(values = c("N ~ 1k" = 0.55, "N ~ 10k+" = 1.0),
                       name = "Population size") +
    labs(x = NULL, y = NULL,
         title = "Neutral simulations - parsimony-saturation summary stats",
         subtitle = "Rows = statistic. Cols = death rate d. Bars = timing x N. Dashed = perfect (=1). Deterministic bars repeat across d columns since deterministic sims have no d parameter.") +
    theme_bw(base_size = 10) +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text.y.left = element_text(angle = 0, face = "bold"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30"))

out4 <- file.path(OUT_DIR, "neutral_summary_grid.png")
ggsave(out4, p_neutral_grid, width = 11, height = 8, dpi = 130)
message("wrote: ", out4)

# ---------------------------------------------------------------------------
# Phylogeny-topology aggregate figures: normalised Robinson-Foulds distance
# between the (subsampled) truth tree and the MEDICC2 reconstruction. RF is
# bounded to [0, 1] by definition (0 = identical topology, 1 = maximally
# different), so both figures share that y-axis.
#
# Saved to a separate directory since they concern phylogeny topology, not
# parsimony-saturation.
# ---------------------------------------------------------------------------
PHY_OUT_DIR <- file.path(ROOT, "figures", "3_inference", "03_MEDICC2",
                         "phylogenies", "aggregate")
dir.create(PHY_OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------ Neutral-only RF grid: one row, cols = d, bars per timing x N --------
neutral_rf <- neutral %>%
    select(timing, N, N_group, d, rf_norm) %>%
    mutate(timing = factor(timing, levels = c("deterministic", "gamma", "markov")))

det_rf <- neutral_rf %>% filter(timing == "deterministic")
non_det_rf <- neutral_rf %>% filter(timing != "deterministic")
det_rf_expanded <- bind_rows(
    det_rf %>% mutate(d = 0.0),
    det_rf %>% mutate(d = 0.5),
    det_rf %>% mutate(d = 0.9)
)
rf_grid_df <- bind_rows(non_det_rf, det_rf_expanded) %>%
    mutate(d_lab = factor(sprintf("d = %.1f", d),
                          levels = c("d = 0.0", "d = 0.5", "d = 0.9")))

p_neutral_rf <- ggplot(rf_grid_df,
                       aes(x = timing, y = rf_norm,
                           fill = timing, alpha = N_group)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             colour = "grey20", linewidth = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    facet_wrap(~ d_lab, nrow = 1) +
    scale_fill_manual(values = timing_colours, name = "Timing") +
    scale_alpha_manual(values = c("N ~ 1k" = 0.55, "N ~ 10k+" = 1.0),
                       name = "Population size") +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "Normalised RF distance",
         title = "Neutral simulations - MEDICC2 vs truth topology (normalised Robinson-Foulds)",
         subtitle = "0 = identical topology (dashed), 1 = maximally different. Bars = timing x N. Deterministic bars repeat across d columns since deterministic sims have no d parameter.") +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30"))

out_rf1 <- file.path(PHY_OUT_DIR, "neutral_rf_grid.png")
ggsave(out_rf1, p_neutral_rf, width = 11, height = 4.5, dpi = 130)
message("wrote: ", out_rf1)

# ------ All-scenario RF: 3 panels (neutral / sel1 / sel2) -------------------
# Reuses the panel()/build_figure() helpers so styling matches the other
# aggregate figures; ref_y = 0 (identical topology) and ylim = [0, 1] since
# RF is bounded by construction.
fig_rf_scen <- build_figure(
    xvar = "truth_EM", yvar = "rf_norm",
    xlab = "Truth E[M] (mean per-cell root->leaf events)",
    ylab = "Normalised RF (MEDICC2 vs truth)",
    big_title = "Tree topology agreement -- normalised Robinson-Foulds distance",
    subtitle  = "0 = identical topology (dashed), 1 = maximally different. Deep truth trees (large E[M]) leave more room for MEDICC2 to disagree.",
    diagonal = FALSE, ref_y = 0, ylim_override = c(0, 1)
)
out_rf2 <- file.path(PHY_OUT_DIR, "rf_by_scenario.png")
ggsave(out_rf2, fig_rf_scen, width = 15, height = 5.5, dpi = 130)
message("wrote: ", out_rf2)
