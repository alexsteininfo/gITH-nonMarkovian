# Plotting constants for the R figures, mirroring plotting_functions.jl.
#
# plotting_functions.jl is the source of truth; this file is kept in sync by hand,
# so a colour changed on one side must be changed on both. Colours are written as
# hex rather than by name because R's colors() has no "teal" -- a name-for-name copy
# would fail silently on the markov panels.

MODEL_COLS <- c(
    deterministic = "#4682B4",   # :steelblue
    gamma         = "#9370DB",   # :mediumpurple
    markov        = "#008080"    # :teal
)

# Used for the sampled lineages overdrawn on the full trees.
HIGHLIGHT_COLS <- c(
    deterministic = "#000080",   # :navy
    gamma         = "#551A8B",   # :purple4
    markov        = "#006400"    # :darkgreen
)

MODEL_LABELS <- c(
    deterministic = "deterministic",
    gamma         = "gamma (k = 5)",
    markov        = "markov"
)

D_LABELS <- c("0" = "d = 0", "0.9" = "d = 0.9")

FS_TITLE <- 16
FS_LABEL <- 14
FS_TICK  <- 12
FS_ANNOT <- 10
