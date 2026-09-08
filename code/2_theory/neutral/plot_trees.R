### Lineage tree figures -- remake of Ycart (2013) Figure 1 from simulated data.
###
### Four 2x3 grids: columns are the timing models, rows are d = 0 / 0.9. The four
### figures form a 2x2 in (tips drawn) x (population sampled from); see
### docs/superpowers/specs/2026-09-08-lineage-tree-figures-design.md.
###
### Reads data/newick/neutral/panels.csv, written by
### code/2_theory/neutral/export_trees_newick.jl. Run that first:
###
###     julia --project=. code/2_theory/neutral/export_trees_newick.jl
###     Rscript code/2_theory/neutral/plot_trees.R

suppressPackageStartupMessages({
    library(ape)
    library(ggtree)
    library(ggplot2)
    library(dplyr)
})

script_path <- normalizePath(
    sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
)
ROOT <- dirname(dirname(dirname(dirname(script_path))))   # code/2_theory/neutral/x.R
source(file.path(ROOT, "code", "paths.R"))
source(file.path(HELPERS, "plotting_functions.R"))

NWK_DIR <- file.path(DATA, "newick", "neutral")
OUT_DIR <- file.path(FIGURES, "2_theory", "neutral")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

### Layout ---------------------------------------------------------------------

# Rectangular dendrogram coordinates for one tree, transposed so that time runs down
# the plot's y axis and tip order across its x axis, both normalised to [0, 1].
#
# ggtree's node x is the cumulative branch length from the root. Under the writer's
# convention (the edge into a node is that cell's lifetime) that makes node x equal to
# the cell's division time, and the founder cell's own lifetime lives in tree$root.edge
# and has to be added back -- otherwise every panel is short by one generation and the
# tips do not reach 1. Same root_off handling as tree_to_plot_data() in
# CNA-PhyloAnalysis/code/2_phylogenies/base_trees.R, which this function follows.
tree_segments <- function(tree, t_end) {
    # suppressWarnings: ggtree passes as.Date / yscale_mapping / hang down to
    # ggplot2::fortify, which warns about unused arguments on ggplot2 >= 3.5. It is a
    # ggtree-internal version mismatch, not something this call can avoid.
    gt       <- suppressWarnings(suppressMessages(ggtree(tree, layout = "rectangular")))$data
    root_off <- if (is.null(tree$root.edge)) 0 else tree$root.edge
    ntip     <- ape::Ntip(tree)

    # ggtree accumulates ~100 parsed branch lengths per tip, so the deepest tips land
    # a few ulps past t_end. The tree is ultrametric by construction (asserted below
    # and in verify_tree_export.jl), so clamp that epsilon -- but only that much: a
    # real overshoot means a wrong t_end and must not be silently squashed into range,
    # where the scale limits would otherwise drop whole branches as "outside the
    # scale range".
    node_t <- (gt$x + root_off) / t_end
    stopifnot(max(node_t) <= 1 + 1e-9)
    node_t <- setNames(pmin(node_t, 1), gt$node)             # -> plot y
    node_i <- setNames((gt$y - 0.5) / ntip, gt$node)         # -> plot x

    edges  <- tree$edge                                      # [parent, child]
    parent <- as.character(edges[, 1])
    child  <- as.character(edges[, 2])

    # a cell's lifetime: vertical, at the child's tip-order position
    life <- data.frame(
        x    = unname(node_i[child]), xend = unname(node_i[child]),
        y    = unname(node_t[parent]), yend = unname(node_t[child]),
        node = edges[, 2]
    )

    # the division: horizontal, at the parent's division time, across its children
    parents <- unique(edges[, 1])
    kids_i  <- split(node_i[child], edges[, 1])
    div <- data.frame(
        x    = unname(vapply(kids_i[as.character(parents)], min, numeric(1))),
        xend = unname(vapply(kids_i[as.character(parents)], max, numeric(1))),
        y    = unname(node_t[as.character(parents)]),
        yend = unname(node_t[as.character(parents)]),
        node = parents
    )

    # the founder cell's own lifetime, above the first division
    root_node <- setdiff(edges[, 1], edges[, 2])
    trunk <- data.frame(
        x    = unname(node_i[as.character(root_node)]),
        xend = unname(node_i[as.character(root_node)]),
        y    = 0, yend = root_off / t_end,
        node = root_node
    )

    rbind(life, div, trunk)
}

# Node ids on the path from the root to any tip whose label is in `ids`.
ancestral_nodes <- function(tree, ids) {
    parent_of <- setNames(tree$edge[, 1], tree$edge[, 2])
    seen      <- logical(max(tree$edge))
    for (tp in which(tree$tip.label %in% ids)) {
        n <- tp
        while (!is.na(parent_of[as.character(n)])) {
            if (seen[n]) break
            seen[n] <- TRUE
            n <- unname(parent_of[as.character(n)])
        }
        seen[n] <- TRUE   # the root itself
    }
    which(seen)
}

### Build one figure -----------------------------------------------------------

panels <- read.csv(file.path(NWK_DIR, "panels.csv"), stringsAsFactors = FALSE)

as_factors <- function(df) {
    df$model_label <- factor(MODEL_LABELS[df$model], levels = unname(MODEL_LABELS))
    df$d_label     <- factor(D_LABELS[df$d],         levels = unname(D_LABELS))
    df
}

build_figure <- function(fig_key) {
    rows <- panels[panels$figure == fig_key, ]
    stopifnot(nrow(rows) == 6)

    segs <- list()
    hits <- list()
    for (i in seq_len(nrow(rows))) {
        r    <- rows[i, ]
        tree <- ape::read.tree(file.path(NWK_DIR, r$newick_file))

        # A malformed export must fail here rather than produce a subtly wrong
        # picture. is.ultrametric is an independent parse of the written file: every
        # leaf is alive at t_end, so every root-to-tip path must have equal length.
        stopifnot(ape::Ntip(tree) == r$n_tips)
        stopifnot(ape::is.ultrametric(tree, tol = 1e-6))

        s       <- tree_segments(tree, r$t_end)
        s$model <- r$model
        s$d     <- as.character(r$d)
        segs[[length(segs) + 1]] <- s

        if (nzchar(r$highlight_file)) {
            ids        <- readLines(file.path(NWK_DIR, r$highlight_file))
            keep       <- s[s$node %in% ancestral_nodes(tree, ids), ]
            hits[[length(hits) + 1]] <- keep
        }
    }
    segs <- as_factors(bind_rows(segs))
    hits <- if (length(hits)) as_factors(bind_rows(hits)) else NULL

    lw <- if (max(rows$n_tips) > 500) 0.15 else 0.4

    ann      <- as_factors(transform(rows, d = as.character(d)))
    ann$text <- sprintf("t_end = %.1f\nn = %d of %d\nmedian depth = %d",
                        rows$t_end, rows$n_tips, rows$N_pop, rows$median_leaf_depth)
    ann$text <- ifelse(rows$repeated == "true",
                       paste0(ann$text, "\n(d = 0 repeated: no death runs)"),
                       ann$text)

    cols <- c(MODEL_COLS,
              setNames(HIGHLIGHT_COLS, paste0(names(HIGHLIGHT_COLS), "_hl")))

    # Where sampled lineages are overdrawn, the whole tree is faded so the sample
    # reads as a distinct layer; a 10% sample keeps nearly the entire backbone, so at
    # full opacity the two layers merge into an unreadable gradient. Figures with no
    # overdraw keep the model colour at full strength.
    base_alpha <- if (is.null(hits)) 1 else 0.3

    p <- ggplot() +
        geom_segment(data = segs,
                     aes(x = x, xend = xend, y = y, yend = yend, colour = model),
                     linewidth = lw, alpha = base_alpha)

    if (!is.null(hits)) {
        p <- p + geom_segment(data = hits,
                              aes(x = x, xend = xend, y = y, yend = yend,
                                  colour = paste0(model, "_hl")),
                              linewidth = lw * 2)
    }

    # vjust is in device space, and scale_y_reverse puts y = 0.02 near the *top* of
    # the panel, so the label has to hang downward from it -- vjust = 0 sends it off
    # the panel and clips all but its last line.
    p +
        geom_text(data = ann, aes(x = 0.02, y = 0.02, label = text),
                  hjust = 0, vjust = 1, size = FS_ANNOT / .pt, inherit.aes = FALSE) +
        facet_grid(d_label ~ model_label) +
        scale_colour_manual(values = cols, guide = "none") +
        scale_y_reverse(limits = c(1.0, 0.0), expand = expansion(mult = 0.03)) +
        scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0.01)) +
        labs(x = NULL, y = expression(t / t[end])) +
        theme_bw(base_size = FS_LABEL) +
        theme(
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid   = element_blank(),
            strip.text   = element_text(size = FS_TITLE),
            axis.text.y  = element_text(size = FS_TICK)
        )
}

### Write ----------------------------------------------------------------------

for (fig_key in unique(panels$figure)) {
    out <- file.path(OUT_DIR, sprintf("trees_%s.png", fig_key))
    ggsave(out, plot = build_figure(fig_key), width = 12, height = 8, dpi = 300)
    cat("->", out, "\n")
}
