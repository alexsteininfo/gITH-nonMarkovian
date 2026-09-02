# Shared machinery for selection scenario 1. Lives here rather than in the three
# `analysis/sim_runs/selection_1/growth_sel1_*.jl` scripts so the accept/retry loop
# exists once instead of three times; those scripts supply only the model-specific
# birth/death distributions and the output path.
#
# The calling script must already have done `using MutationLoadDynamics,
# Distributions, Random, Statistics` and included `types_selection1.jl`.

"""
    ncritic_grid(N_target) -> Vector{Int}

Ten log-spaced injection sizes from 1 to `N_target ÷ 2`.

Log spacing makes the injection *times* approximately equally spaced: during
exponential growth `N(t) ≈ e^{rt}`, so `t ≈ ln(N)/r`. Ten values rather than
twenty because twenty log-spaced integers over `[1, 500]` collide at the low end
(1, 1, 2, 3, 4, 5, …) and dedup down to about seventeen.
"""
ncritic_grid(N_target::Int) =
    unique(round.(Int, exp.(range(0.0, log(N_target ÷ 2), length = 10))))
