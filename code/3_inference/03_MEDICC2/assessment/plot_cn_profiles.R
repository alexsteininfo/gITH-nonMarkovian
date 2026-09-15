# plot_cn_profiles.R
#
# (A) Allele-specific CN profiles for a stratified + random sample of 6 cells.
# For single cells the MEDICC2 leaf profile equals its input (truth), so we
# only draw the MEDICC2 profile; alleles A/B are overlaid on the same panel
# with a small +/- y-offset (and different colours) for readability.
#
# Inputs:
#   data/CN_subsampled/<param_dir>/sim1_tree.nwk
#   data/CN_subsampled/<param_dir>/sim1_truth_events.tsv   (for burden-based sampling only)
#   data/MEDICC2/treeinference/<sim_rel>/sim1_final_cn_profiles.tsv
#
# Outputs:
#   figures/3_inference/03_MEDICC2/assessment/cn_profiles__<sim_slug>.png
#
# Usage:
#   Rscript code/3_inference/03_MEDICC2/assessment/plot_cn_profiles.R \
#       [--sim <sim_rel>] [--outdir <abs_path_or_repo_rel>]
#
# Default sim_rel: neutral/gamma/neutral_gamma_N1000_d0.5_k5.0_n100/sim1
# Default outdir : <ROOT>/figures/3_inference/03_MEDICC2/assessment/

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(readr)
    library(ggplot2); library(cowplot)
})

args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE))
if (length(this_file) != 1L) stop("cannot resolve --file")
ROOT <- normalizePath(file.path(dirname(this_file), "..", "..", "..", ".."))

source(file.path(ROOT, "code", "X_helpers", "medicc2_io.R"))
source(file.path(ROOT, "code", "X_helpers", "plotting_functions.R"))

# --- Arg parsing (plain commandArgs, no optparse) ---
user_args <- commandArgs(trailingOnly = TRUE)
sim_rel <- "neutral/gamma/neutral_gamma_N1000_d0.5_k5.0_n100/sim1"
outdir_arg <- NULL
i <- 1
while (i <= length(user_args)) {
    if (user_args[i] == "--sim") {
        sim_rel <- user_args[i + 1L]; i <- i + 2L
    } else if (user_args[i] == "--outdir") {
        outdir_arg <- user_args[i + 1L]; i <- i + 2L
    } else {
        stop("unknown arg: ", user_args[i])
    }
}
message("sim: ", sim_rel)

TRUTH_DIR   <- file.path(ROOT, "data", "CN_subsampled", dirname(sim_rel))
MEDICC2_DIR <- file.path(ROOT, "data", "MEDICC2", "treeinference", sim_rel)
# Absolute --outdir wins; a relative --outdir is resolved against ROOT; default
# preserves the pre-existing behaviour (assessment/ next to the other plots).
OUT_DIR <- if (is.null(outdir_arg)) {
    file.path(ROOT, "figures", "3_inference", "03_MEDICC2", "assessment")
} else if (startsWith(outdir_arg, "/")) {
    outdir_arg
} else {
    file.path(ROOT, outdir_arg)
}
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Load ---
# truth_profiles is not needed: for a leaf cell, MEDICC2's reconstructed
# profile is identical to the truth (truth was its input), so we only need
# the MEDICC2 per-bin profiles. Truth tree + events are still needed to
# compute per-cell burden for the stratified cell sample.
truth_tree_raw <- read_truth_tree(TRUTH_DIR)
truth_events   <- read_truth_events(TRUTH_DIR)
med_profiles   <- read_medicc2_profiles(MEDICC2_DIR)

message("truth cells: ", length(truth_tree_raw$tip.label))
message("medicc2 unique sample_ids: ", n_distinct(med_profiles$sample_id))

# --- Pick 4 stratified + 2 random cells by truth burden ---
collapsed <- collapse_unary(truth_tree_raw)
bt <- per_cell_burden_truth(collapsed, truth_events)
picks <- stratified_cell_sample(bt, n_stratified = 4L, n_random = 2L, seed = 42L)
message("selected cells: ", paste(picks, collapse = ", "))

# Facet order = ascending truth burden (fills row-wise across 2 cols).
cell_order <- bt %>% filter(cell %in% picks) %>% arrange(M) %>% pull(cell)

# --- MEDICC2 per-bin profiles: keep picked cells, drop chrX ---
med_long <- med_profiles %>%
    filter(sample_id %in% picks, chrom != "chrX") %>%
    select(cell = sample_id, chrom, start, stop = end, cn_a, cn_b) %>%
    pivot_longer(c(cn_a, cn_b), names_to = "allele", values_to = "cn")

message("medicc2 rows after filter: ", nrow(med_long))

# --- Compute a single per-cell x-coordinate spanning all chromosomes ---
# Chromosomes are numeric 1..22 after dropping chrX; order them numerically.
chroms_in_data <- unique(med_long$chrom)
chr_num <- suppressWarnings(as.integer(sub("^chr", "", chroms_in_data)))
chr_order <- chroms_in_data[order(chr_num, chroms_in_data)]

chr_span <- med_long %>%
    group_by(chrom) %>% summarise(len = max(stop), .groups = "drop") %>%
    mutate(chrom = factor(chrom, levels = chr_order)) %>%
    arrange(chrom) %>%
    mutate(offset = cumsum(as.numeric(len)) - as.numeric(len))

offset_lookup <- setNames(chr_span$offset, as.character(chr_span$chrom))

med_plot <- med_long %>%
    mutate(chrom = factor(chrom, levels = chr_order),
           x_start = as.numeric(start) + offset_lookup[as.character(chrom)],
           x_stop  = as.numeric(stop)  + offset_lookup[as.character(chrom)])

# Chromosome divider positions and label midpoints
chr_dividers <- chr_span$offset[-1]
chr_labels <- chr_span %>%
    mutate(mid = offset + as.numeric(len) / 2) %>%
    select(chrom, mid)

# --- Overlay alleles on the same panel with a small +/- y-offset ---
# Single-cell CN values are integers, so a shift of 0.05 keeps 0/1/2/3 clearly
# separable while making the two allele tracks visible side by side.
y_offset <- 0.05
med_plot <- med_plot %>%
    mutate(y = cn + if_else(allele == "cn_a", y_offset, -y_offset),
           cell = factor(cell, levels = cell_order),
           allele = factor(allele, levels = c("cn_a", "cn_b"),
                           labels = c("Allele A", "Allele B")))

allele_colours <- c("Allele A" = "#1F77B4", "Allele B" = "#D62728")

p <- ggplot(med_plot) +
    geom_segment(aes(x = x_start, xend = x_stop, y = y, yend = y, colour = allele),
                 linewidth = 0.9) +
    geom_vline(xintercept = chr_dividers, colour = "grey85", linewidth = 0.3) +
    facet_wrap(vars(cell), ncol = 2L) +
    scale_x_continuous(breaks = chr_labels$mid,
                       labels = sub("^chr", "", chr_labels$chrom),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = 0:4) +
    scale_colour_manual(values = allele_colours, name = NULL) +
    labs(x = "Chromosome", y = "Copy number",
         title = "MEDICC2 CN profiles (Allele A: +0.05, Allele B: -0.05 offset)",
         subtitle = sim_rel) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.spacing.x = unit(0.4, "lines"),
          panel.spacing.y = unit(0.4, "lines"),
          strip.text = element_text(face = "bold"),
          legend.position = "top")

out_path <- file.path(OUT_DIR, paste0("cn_profiles__", sim_slug(sim_rel), ".png"))
ggsave(out_path, p, width = 12, height = 9, dpi = 120)
message("wrote: ", out_path)
