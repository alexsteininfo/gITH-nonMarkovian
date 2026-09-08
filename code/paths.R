# Project paths for R scripts, mirroring code/paths.jl.
#
# The calling script defines ROOT itself -- by climbing dirname() from its own file
# path -- and then sources this. Same division of labour as the Julia side, and for
# the reason CLAUDE.md gives: a helper that guessed the depth on its caller's behalf
# is exactly how the wrong environment got activated in the previous layout.
#
#     source(file.path(ROOT, "code", "paths.R"))

if (!exists("ROOT")) {
    stop("paths.R: ROOT must be defined before sourcing this file")
}
if (!file.exists(file.path(ROOT, "Project.toml"))) {
    stop(sprintf("ROOT = %s has no Project.toml -- was this script moved?", ROOT))
}

DATA    <- file.path(ROOT, "data")
FIGURES <- file.path(ROOT, "figures")
HELPERS <- file.path(ROOT, "code", "X_helpers")
