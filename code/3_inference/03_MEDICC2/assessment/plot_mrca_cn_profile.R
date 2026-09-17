# plot_mrca_cn_profile.R
#
# Allele-specific CN profiles at the sample-MRCA for one sim: truth as a thick
# semi-transparent underlay, MEDICC2 as a thin opaque overlay. Same chromosome
# layout / allele +/- y-offset convention as plot_cn_profiles.R.
#
# The truth is drawn from its native variable-length segments (one line per
# maximal same-CN run per haplotype); MEDICC2 is drawn from its 1 Mb bin grid.
# Any spurious CNA at the MEDICC2 MRCA shows up as an opaque bin where the
# truth line stays flat at diploid.
#
# Default sim = neutral/markov/neutral_markov_N1000_d0.9_n100/sim1 -- the case
# with 27 hallucinated chr21-arm gain bins at the MRCA. Three other sims in
# the current sweep also have imperfect MRCA reconstructions (see
# data/MEDICC2/benchmark/mrca_reconstruction.csv):
#   neutral/markov/neutral_markov_N10000_d0.9_n100/sim1      (1 halluc hbin)
#   selection_1/gamma/sel1_gamma_N1000_d0.9_k5.0_s0.1_n100/sim1 (3 halluc hbins)
#   selection_1/gamma/sel1_gamma_N10000_d0.5_k5.0_s1.7_n100/sim1 (6 missed hbins)
#
# Usage:
#   Rscript code/3_inference/03_MEDICC2/assessment/plot_mrca_cn_profile.R \
#       [--sim <sim_rel>] [--outdir <abs_path_or_repo_rel>]

suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(readr); library(ggplot2); library(ape)
})

args_all <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args_all, value = TRUE))
if (length(this_file) != 1L) stop("cannot resolve --file")
ROOT <- normalizePath(file.path(dirname(this_file), "..", "..", "..", ".."))

source(file.path(ROOT, "code", "X_helpers", "medicc2_io.R"))

# --- Args ------------------------------------------------------------------
user_args <- commandArgs(trailingOnly = TRUE)
sim_rel <- "neutral/markov/neutral_markov_N1000_d0.9_n100/sim1"
outdir_arg <- NULL
i <- 1L
while (i <= length(user_args)) {
    if      (user_args[i] == "--sim")    { sim_rel   <- user_args[i + 1L]; i <- i + 2L }
    else if (user_args[i] == "--outdir") { outdir_arg <- user_args[i + 1L]; i <- i + 2L }
    else stop("unknown arg: ", user_args[i])
}
message("sim: ", sim_rel)

TRUTH_DIR   <- file.path(ROOT, "data", "CN_subsampled", dirname(sim_rel))
MEDICC2_DIR <- file.path(ROOT, "data", "MEDICC2", "treeinference", sim_rel)
OUT_DIR <- if (is.null(outdir_arg)) {
    file.path(ROOT, "figures", "3_inference", "03_MEDICC2", "MRCA_test")
} else if (startsWith(outdir_arg, "/")) {
    outdir_arg
} else {
    file.path(ROOT, outdir_arg)
}
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Load MRCA labels via the same collapse/drop-diploid convention -------
truth_raw <- read_truth_tree(TRUTH_DIR)
med_full  <- read_medicc2_tree(MEDICC2_DIR)
med_tree  <- if ("diploid" %in% med_full$tip.label)
    ape::drop.tip(med_full, "diploid") else med_full

truth_tree <- collapse_unary(truth_raw)$tree
truth_mrca_name <- truth_tree$node.label[1L]
med_mrca_name   <- med_tree$node.label[1L]
message(sprintf("MRCA labels: truth=%s, medicc=%s", truth_mrca_name, med_mrca_name))

# --- CN profiles at the MRCA (chrX already dropped by the helper) ---------
truth_prof <- read_truth_profiles(TRUTH_DIR) %>%
    filter(name == truth_mrca_name) %>%
    select(chrom, haplotype, start, stop, cn) %>%
    mutate(allele = ifelse(haplotype == 1L, "Allele A", "Allele B"))

med_prof <- read_medicc2_profiles(MEDICC2_DIR) %>%
    filter(sample_id == med_mrca_name, chrom != "chrX") %>%
    select(chrom, start, stop = end, cn_a, cn_b) %>%
    pivot_longer(c(cn_a, cn_b), names_to = "allele_key", values_to = "cn") %>%
    mutate(allele = ifelse(allele_key == "cn_a", "Allele A", "Allele B"))

# --- Chromosome axis: same conventions as plot_cn_profiles.R -------------
chroms_in_data <- unique(med_prof$chrom)
chr_num <- suppressWarnings(as.integer(sub("^chr", "", chroms_in_data)))
chr_order <- chroms_in_data[order(chr_num, chroms_in_data)]

chr_span <- med_prof %>%
    group_by(chrom) %>% summarise(len = max(stop), .groups = "drop") %>%
    mutate(chrom = factor(chrom, levels = chr_order)) %>%
    arrange(chrom) %>%
    mutate(offset = cumsum(as.numeric(len)) - as.numeric(len))

offset_lookup <- setNames(chr_span$offset, as.character(chr_span$chrom))

add_x_coords <- function(df) {
    df %>%
        filter(chrom %in% chr_order) %>%
        mutate(chrom = factor(chrom, levels = chr_order),
               x_start = as.numeric(start) + offset_lookup[as.character(chrom)],
               x_stop  = as.numeric(stop)  + offset_lookup[as.character(chrom)])
}

y_offset <- 0.05
truth_plot <- add_x_coords(truth_prof) %>%
    mutate(y = cn + if_else(allele == "Allele A", y_offset, -y_offset),
           allele = factor(allele, levels = c("Allele A", "Allele B")))

med_plot <- add_x_coords(med_prof) %>%
    mutate(y = cn + if_else(allele == "Allele A", y_offset, -y_offset),
           allele = factor(allele, levels = c("Allele A", "Allele B")))

chr_dividers <- chr_span$offset[-1]
chr_labels <- chr_span %>%
    mutate(mid = offset + as.numeric(len) / 2) %>%
    select(chrom, mid)

allele_colours <- c("Allele A" = "#1F77B4", "Allele B" = "#D62728")

# --- Plot: truth thick + transparent underlay, MEDICC2 thin opaque overlay
p <- ggplot() +
    geom_vline(xintercept = chr_dividers, colour = "grey85", linewidth = 0.3) +
    geom_segment(data = truth_plot,
                 aes(x = x_start, xend = x_stop, y = y, yend = y, colour = allele),
                 linewidth = 3.2, alpha = 0.28,
                 lineend = "butt") +
    geom_segment(data = med_plot,
                 aes(x = x_start, xend = x_stop, y = y, yend = y, colour = allele),
                 linewidth = 0.9, alpha = 1.0,
                 lineend = "butt") +
    scale_x_continuous(breaks = chr_labels$mid,
                       labels = sub("^chr", "", chr_labels$chrom),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = 0:4) +
    scale_colour_manual(values = allele_colours, name = NULL) +
    labs(x = "Chromosome", y = "Copy number",
         title = sprintf("Sample-MRCA CN profile: truth (thick, transparent) vs MEDICC2 (thin, opaque)"),
         subtitle = sprintf("%s   |   truth MRCA = %s, MEDICC2 MRCA = %s   |   Allele A: +0.05, Allele B: -0.05 offset",
                            sim_rel, truth_mrca_name, med_mrca_name)) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "top",
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(colour = "grey30", size = 9))

out_path <- file.path(OUT_DIR, paste0("mrca_cn_profile__", sim_slug(sim_rel), ".png"))
ggsave(out_path, p, width = 12, height = 4.2, dpi = 130)
message("wrote: ", out_path)
