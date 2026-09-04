# Sample-size table, seed derivation, atomic writes and the shard/file drivers for
# subsampling. The draw itself is `MutationLoadDynamics.sample_leaves`.
#
# Requires, from the caller: `using MutationLoadDynamics` and `SubsampleResult`
# (helpers/types_subsampled.jl) plus the scenario's params types.

using AbstractTrees
using Random
using Serialization

# ── Atomic writes ─────────────────────────────────────────────────────────────

"""
    serialize_atomic(path, data)

Serialize `data` to `path` via a temporary file plus rename, instead of
`serialize(path, data)` directly.

Every stage in this pipeline (this file's `subsample_file` and the three
`process_*_subsampled.jl` scripts) decides "already done" purely from
`isfile(path)`, so a resumed run treats *any* file at `path` as complete.
`serialize` writes in place, so a kill mid-write — a real hazard on a
multi-hour sweep — would leave a truncated-but-present file that a later run
silently skips instead of regenerating. Writing to `path * ".tmp"` and then
renaming into place means `path` only ever exists once the write is complete;
`mv` within the same directory is atomic on the local filesystem.
"""
function serialize_atomic(path::AbstractString, data)
    tmp = path * ".tmp"
    serialize(tmp, data)
    mv(tmp, path; force = true)
    return nothing
end

# ── Sample sizes ──────────────────────────────────────────────────────────────

"""
    sample_sizes(N_target) -> Vector{Int}

Cell counts to draw from a tree grown to `N_target`, largest first.

An explicit table rather than `round.(fracs .* N_target)`, so the values in the
data can never drift with a change of rounding convention, and so a newly
introduced `N_target` fails loudly instead of silently acquiring sizes nobody
chose. `102 = round(0.1 * 1024)`, `1638 = round(0.1 * 16384)`,
`164 = round(0.01 * 16384)`.
"""
function sample_sizes(N_target::Int)
    N_target == 1_000  && return [100]
    N_target == 1_024  && return [102]
    N_target == 10_000 && return [1_000, 100]
    N_target == 16_384 && return [1_638, 164]
    error("no sample-size rule for N_target = $N_target — add one to sample_sizes")
end

# ── Seeding ───────────────────────────────────────────────────────────────────

"""
    sample_seed(stem, sim_index, n) -> UInt64

Seed for one draw. Depends only on values recorded in the output file, so any
single draw replays in isolation, and each threaded task builds its own rng, so
threading cannot perturb the samples.
"""
sample_seed(stem::AbstractString, sim_index::Int, n::Int) = hash((stem, sim_index, n))

# ── Shard driver ──────────────────────────────────────────────────────────────

"""
    subsample_shard(sims, stem, n) -> Vector{SubsampleResult{P}}

Draw one `n`-cell subsample from every simulation in one raw shard.

`P` is read off the source result type, so the same function serves all three
scenarios: `fieldtype(GrowthSimResult, :params)` is `SimParams`, and likewise for
`Sel1SimResult` and `Sel2SimResult`. The returned vector is concretely typed, so the
serialized file deserializes as a `Vector{SubsampleResult{P}}`.

Simulations whose `tree_root` is `nothing` are dropped and source order is kept,
matching `analysis/processing/process_*.jl` exactly — so index `i` of the output
lines up with index `i` of the corresponding array under `data/processed/`.
"""
function subsample_shard(sims::Vector{T}, stem::AbstractString, n::Int) where {T}
    P = fieldtype(T, :params)
    R = SubsampleResult{P}

    out = Vector{Union{Nothing, R}}(nothing, length(sims))
    Threads.@threads for i in eachindex(sims)
        sim = sims[i]
        isnothing(sim.tree_root) && continue
        seed = sample_seed(stem, i, n)
        # The draw itself lives in MutationLoadDynamics.jl (v0.3.0+). Only the
        # record type, the seed derivation and the file layout are ours — see
        # docs/superpowers/specs/2026-09-04-repo-architecture-design.md. Do not
        # move `SubsampleResult` into the package: it is defined in `Main`, and
        # 525 serialized shards resolve it by that module path.
        s = sample_leaves(sim.tree_root, n; seed = seed)
        out[i] = R(s.root, sim.params, n, s.N_full, i, seed, s.sampled_ids)
    end
    return R[r for r in out if !isnothing(r)]
end

"""
    subsample_file(raw_path, out_dir, ns)

Deserialize one raw shard and write one subsampled shard per entry of `ns`, named
`<stem>_n<n>.jls` under `out_dir`.

The shard is read once for all sample sizes — deserializing the trees dominates the
runtime. Resumable in the style of `analysis/processing/process_sel1.jl`: a sample
size whose output already exists is skipped, so an interrupted sweep resumes and a
forced redo is a file deletion. Because the seed is a pure function of
`(stem, sim_index, n)`, a redo reproduces the original draws exactly.
"""
function subsample_file(raw_path::String, out_dir::String, ns::Vector{Int})
    stem = splitext(basename(raw_path))[1]
    todo = [n for n in ns if !isfile(joinpath(out_dir, "$(stem)_n$(n).jls"))]

    if isempty(todo)
        println("  Skipping (already subsampled): $(basename(raw_path))")
        return
    end

    println("  Subsampling: $(basename(raw_path))  n ∈ $todo")
    sims = deserialize(raw_path)
    mkpath(out_dir)

    for n in todo
        res     = subsample_shard(sims, stem, n)
        skipped = length(sims) - length(res)
        skipped > 0 &&
            @warn "  $stem: skipped $skipped / $(length(sims)) sims with nothing tree"
        # Atomic: the `isfile` resumability check above must never see a truncated file.
        serialize_atomic(joinpath(out_dir, "$(stem)_n$(n).jls"), res)
        println("    → n=$n: $(length(res)) subsampled trees")
    end

    # Raw shards run to ~400 MB; drop the trees before moving to the next file.
    sims = nothing
    GC.gc()
    return
end
