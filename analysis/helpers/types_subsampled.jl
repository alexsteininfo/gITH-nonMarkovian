# Data structure for subsampled lineage trees.
#
# A `SubsampleResult` records one random draw of `n` cells from one simulated tree,
# together with the induced lineage tree: the sampled leaves plus *every* ancestor
# of a sampled leaf, with unary nodes retained rather than collapsed. That is the
# same shape `prune_tree!` leaves behind when a lineage dies out, so the stored root
# is a `BinaryNode{NonMarkovCell}` indistinguishable in kind from a full tree and
# every helper in `tree_analysis.jl` applies to it unchanged.
#
# `include`d (not imported) by every subsampling stage so that independent scripts
# see identical definitions. Deliberately defines *only* `SubsampleResult`: each
# consuming script also includes the scenario params types it needs
# (`types.jl` / `types_selection1.jl` / `types_selection2.jl`), following the
# convention of `analysis/processing/process_sel1.jl`. BinaryNode and NonMarkovCell
# are exported from MutationLoadDynamics.

"""
One random subsample of one simulation.

`params` is the source simulation's parameter struct verbatim — `SimParams`,
`Sel1Params` or `Sel2Params` — so a subsampled shard carries its own design point.

`trajectory` is not stored: it is a population-level time series, unaffected by
sampling, and already in `data/raw/`. Scenario-specific records (`Sel1Injection`,
`n_restarts`) are not stored either, for the same reason; the processing stage
copies them across from `data/processed/`.

# Fields
- `tree_root` — induced tree: sampled leaves plus all their ancestors
- `params` — the source simulation's parameters
- `n` — cells sampled
- `N_full` — leaf count of the source tree
- `sim_index` — index into the *unfiltered* source `Vector`, so a subsample can be
  traced back to its raw simulation even though `nothing`-tree entries are dropped
- `seed` — rng seed of this draw; the draw is a pure function of
  `(tree, seed, n)` and so replays in isolation
- `sampled_ids` — `NonMarkovCell.id` of each drawn cell, in draw order

Forward-compatibility note: there is no draw-replicate field. `sample_seed`
currently makes `(stem, sim_index, n)` determine the draw uniquely, so "k
independent draws per (sim, n)" cannot be added later by widening this struct —
doing so would make the 525 existing shards unreadable under `Serialization`
(see `CLAUDE.md`). It would have to come either as a new struct/filename
convention, or by overloading `n` in the filename (e.g. one draw's `n` tagged
with a replicate suffix) rather than as a field here.
"""
struct SubsampleResult{P}
    tree_root::BinaryNode{NonMarkovCell}
    params::P
    n::Int
    N_full::Int
    sim_index::Int
    seed::UInt64
    sampled_ids::Vector{Int64}
end
