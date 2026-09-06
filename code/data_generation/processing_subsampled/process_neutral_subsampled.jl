using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using AbstractTrees
using Serialization
using Random

include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "subsampling.jl"))   # for the sample_sizes table

# ── Paths ─────────────────────────────────────────────────────────────────────

const RAW  = joinpath(DATA, "raw_subsampled",       "neutral")
const PROC = joinpath(DATA, "processed_subsampled", "neutral")

# ── Processing ────────────────────────────────────────────────────────────────
#
# Stage 2b: the same quantities `process_neutral.jl` extracts from a full tree,
# recomputed on a random subsample of n cells.
#
#   params       — SimParams of the source simulation (nu = 2.0, s = 0.0)
#   mut_per_cell — per-cell burden of each sampled cell (Vector{Int}, length n)
#   sfs          — sfs[k] = mutations carried by exactly k of the n sampled cells
#   leaf_depths  — divisions from founder to each sampled cell
#
# `mut_per_cell` and `leaf_depths` are properties of individual cells and are
# *identical* to the values in `data/processed/` for those same cells: the induced
# tree keeps every ancestor, so each cell's path to the founder is unchanged. They
# are stored because a subsample of a distribution is not the distribution, and
# because their invariance is what `verify_subsampling.jl` checks. `sfs` is the
# quantity sampling genuinely distorts, and the point of the exercise.
#
# `sfs` has length n, not N — the deepest bin `sfs[n]` holds the mutations carried by
# every sampled cell, including the founder's own.
#
# Order is inherited from the subsampled shard, which inherited it from the raw
# shard under the same nothing-tree filter as `process_neutral.jl`. So index i of
# every array here is the same simulation as index i under `data/processed/neutral/`.

const SUBDIRS = ("params", "mut_per_cell", "sfs", "leaf_depths")

function process_file(raw_path::String, proc_dir::String, label::String)
    stem = splitext(basename(raw_path))[1]

    # Resumable, like the subsampling stage: delete the processed file to force a redo.
    if all(isfile(joinpath(proc_dir, sub, stem * ".jls")) for sub in SUBDIRS)
        println("  Skipping (already processed): $(basename(raw_path))")
        return
    end

    println("  Processing: $(basename(raw_path))")

    subs = deserialize(raw_path)::Vector{SubsampleResult{SimParams}}

    all_params       = SimParams[]
    all_mut_per_cell = Vector{Int}[]
    all_sfs          = Vector{Int}[]
    all_leaf_depths  = Vector{Int}[]

    for s in subs
        root   = s.tree_root
        depths = leaf_depths(root)

        # The stored n is authoritative for the sfs length; a mismatch means the
        # shard and its filename have come apart, which must not be averaged over.
        length(depths) == s.n ||
            error("$label: sim $(s.sim_index) has $(length(depths)) leaves, expected n = $(s.n)")

        push!(all_params,       s.params)
        push!(all_mut_per_cell, mutations_per_cell(root))
        push!(all_sfs,          sitefrequencyspectrum(root, s.n))
        push!(all_leaf_depths,  depths)
    end

    for (subdir, data) in (
        ("params",        all_params),
        ("mut_per_cell",  all_mut_per_cell),
        ("sfs",           all_sfs),
        ("leaf_depths",   all_leaf_depths),
    )
        outdir = joinpath(proc_dir, subdir)
        mkpath(outdir)
        # Atomic: the `isfile`-based resumability check above must never see a
        # truncated file left by a kill mid-write.
        serialize_atomic(joinpath(outdir, stem * ".jls"), data)
    end

    println("    → saved $(length(all_params)) results to $proc_dir")
end

# The grids mirror `subsample_neutral.jl`, with `_n$(n)` appended.

println("\n── Gamma model ────────────────────────────────────────────────────────")
for N in (1_000, 10_000), d in (0.0, 0.5, 0.9), n in sample_sizes(N)
    fname    = "neutral_gamma_N$(N)_d$(d)_k5.0_n$(n).jls"
    raw_path = joinpath(RAW, "gamma", fname)
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, joinpath(PROC, "gamma"), "gamma N=$N d=$d n=$n")
end

println("\n── Markov model ───────────────────────────────────────────────────────")
for N in (1_000, 10_000), d in (0.0, 0.5, 0.9), n in sample_sizes(N)
    fname    = "neutral_markov_N$(N)_d$(d)_n$(n).jls"
    raw_path = joinpath(RAW, "markov", fname)
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, joinpath(PROC, "markov"), "markov N=$N d=$d n=$n")
end

println("\n── Deterministic model ────────────────────────────────────────────────")
for N in (1_024, 16_384), n in sample_sizes(N)
    fname    = "neutral_deterministic_N$(N)_n$(n).jls"
    raw_path = joinpath(RAW, "deterministic", fname)
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, joinpath(PROC, "deterministic"), "deterministic N=$N n=$n")
end

println("\nAll done.")
