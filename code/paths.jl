# Project paths, resolved from the ROOT that the including script has already
# computed and checked. No script under code/ counts directory levels to reach
# data or figures — that counting is what silently activated the wrong environment
# in the previous layout (see CLAUDE.md).
#
# Included by every executable script under code/, immediately after its header:
#     include(joinpath(ROOT, "code", "paths.jl"))

const DATA    = joinpath(ROOT, "data")
const FIGURES = joinpath(ROOT, "figures")
const HELPERS = joinpath(ROOT, "code", "X_helpers")
