# Shared machinery for stage 4: copy-number evolution on the subsampled trees.
#
# Requires, from the caller: `using CopyNumberEvolution, MutationLoadDynamics`,
# plus `SubsampleResult` (helpers/types_subsampled.jl) already included.
#
# See docs/superpowers/specs/2026-09-07-cn-evolution-design.md for the model and
# scope rationale.

using Serialization

# ── Model ─────────────────────────────────────────────────────────────────────
#
# hg38 female assembly, 1Mb bins (coarsened from the 500kb DLP+ default to keep
# this first pass's disk footprint small). Rate is FromEdgeMutations() — exact
# identity: the edges of a MutationLoadDynamics-derived tree already carry the
# correct division-level mutation count via
# CopyNumberEvolutionMutationLoadDynamicsExt, so one mutation becomes one CNA with
# no extra randomness. Every other CNAModel component is left at its package
# default: UniformChromosome target, ExtentMixture() (focal-only — no arm- or
# chromosome-level events), GainLoss(0.5), NoWGD(), RejectAndRedraw() viability,
# Diploid() root (clean start, no truncal alterations).

const ASSEMBLY = hg38(:female)
const GRID     = BinGrid(ASSEMBLY, 1_000_000)
const CN_MODEL = CNAModel(rate = FromEdgeMutations())

# ── Scope ─────────────────────────────────────────────────────────────────────
#
# First pass: one simulation per shard, not every sim in it. Raise this once
# runtime/disk at the current scope is understood — see the design doc.
const SIMS_PER_SHARD = 1

# ── Seeding ───────────────────────────────────────────────────────────────────

"""
    cn_seed(stem, sim_index) -> Int64

Seed for one simulation's copy-number evolution. A pure function of values already
recorded in the raw_subsampled filename (`stem`) and the `SubsampleResult` itself
(`sim_index`), so any run reproduces exactly. The hash is masked to fit in Int64.
"""
cn_seed(stem::AbstractString, sim_index::Int) = Int64(hash((stem, sim_index, :cn_evolution)) & typemax(Int64))

# ── Per-simulation driver ────────────────────────────────────────────────────

"""
    cn_evolve_sim(s::SubsampleResult, stem, out_dir) -> Bool

Simulate copy-number alterations on one subsampled tree and write MEDICC2-ready
output plus ground truth into `out_dir`. Returns `false` without doing anything if
every output file for this `sim_index` already exists (resumable), `true` if it ran.
"""
function cn_evolve_sim(s::SubsampleResult, stem::AbstractString, out_dir::String)
    prefix  = joinpath(out_dir, "sim$(s.sim_index)_")
    outputs = (prefix * "cells.tsv", prefix * "truth_profiles.tsv",
               prefix * "truth_events.tsv", prefix * "tree.nwk")
    all(isfile, outputs) && return false

    tree = PhyloTree(s.tree_root)
    seed = cn_seed(stem, s.sim_index)
    res  = simulate_cnas(tree, ASSEMBLY, CN_MODEL; seed = seed)
    mat  = CNMatrix(res, GRID)

    mkpath(out_dir)
    write_medicc2(prefix * "cells.tsv", mat)
    write_profiles(prefix * "truth_profiles.tsv", res)
    write_events(prefix * "truth_events.tsv", res)
    write_newick(prefix * "tree.nwk", tree; branchlength = :divisions)
    return true
end

# ── Per-shard driver ──────────────────────────────────────────────────────────

"""
    cn_evolve_file(raw_path, out_dir)

Deserialize one `data/raw_subsampled/` shard and run copy-number evolution on its
first `SIMS_PER_SHARD` simulations (in stored order — the order `raw_subsampled/`
already filters and preserves from the raw simulations), writing into
`out_dir/<stem>/`.
"""
function cn_evolve_file(raw_path::String, out_dir::String)
    stem    = splitext(basename(raw_path))[1]
    sim_dir = joinpath(out_dir, stem)

    println("  Evolving: $(basename(raw_path))")
    subs = deserialize(raw_path)
    todo = subs[1:min(SIMS_PER_SHARD, length(subs))]

    n_done = 0
    for s in todo
        cn_evolve_sim(s, stem, sim_dir) && (n_done += 1)
    end
    println(n_done == 0 ? "    → already done" : "    → $n_done new simulation(s) evolved")
    return nothing
end
