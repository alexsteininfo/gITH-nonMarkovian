# Convert MutationLoadDynamics.jl simulator outputs into EvoTree, preserving the
# provenance keys the EvoTracerMutationLoadDynamicsExt extension writes.
#
# `include`d (not imported) by the six mutation-rate and six selection drivers,
# so that Serialization sees identical definitions across independent scripts.

using MutationLoadDynamics
using EvoTracer

"""
    evotree_from_full(sim; primary = :mutations, label_prefix = "cell") -> EvoTree

Wrap `EvoTree(sim.tree_root; primary, label_prefix)`. `sim` is any of
`GrowthSimResult`, `Sel1SimResult`, `Sel2SimResult` — all three carry `tree_root`.
"""
function evotree_from_full(sim; primary::Symbol = :mutations,
                           label_prefix::AbstractString = "cell")
    return EvoTree(sim.tree_root; primary = primary, label_prefix = label_prefix)
end

"""
    evotree_from_subsample(sub::SubsampleResult; primary = :mutations,
                           label_prefix = "cell") -> EvoTree

Convert a `SubsampleResult`. The extension only defines
`EvoTree(::LeafSample; …)`; `SubsampleResult` is a repo-local struct with the
same fields split apart, so this rebuilds the same provenance keys manually.
"""
function evotree_from_subsample(sub::SubsampleResult; primary::Symbol = :mutations,
                                label_prefix::AbstractString = "cell")
    t = EvoTree(sub.tree_root; primary = primary, label_prefix = label_prefix)
    p = provenance(t)
    p[:source]       = "MutationLoadDynamics LeafSample (via SubsampleResult)"
    p[:n_sampled]    = sub.n
    p[:n_population] = sub.N_full
    p[:rho]          = sub.N_full > 0 ? sub.n / sub.N_full : NaN
    p[:seed]         = sub.seed
    p[:replicate]    = sub.sim_index
    return t
end

"""
    prepared_full(sim; primary = :mutations) -> EvoTree

Composed with `prepare_tree(_; extract_clade = false)` — the shape every driver
consumes.
"""
prepared_full(sim; primary::Symbol = :mutations) =
    prepare_tree(evotree_from_full(sim; primary = primary); extract_clade = false)

"""
    prepared_sub(sub; primary = :mutations) -> EvoTree
"""
prepared_sub(sub; primary::Symbol = :mutations) =
    prepare_tree(evotree_from_subsample(sub; primary = primary); extract_clade = false)
