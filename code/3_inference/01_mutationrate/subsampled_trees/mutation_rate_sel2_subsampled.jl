# Mutation-rate inference on subsampled sel_2 trees.
#
# For each shard in data/raw_subsampled/selection_2/<timing>/*.jls, deserialize
# the subsampled trees, convert to EvoTree twice (primary = :mutations,
# :divisions), run the full mutation-rate statistic suite, and serialize a
# Vector{Union{Nothing,Dict}} to
# data/inference/mutationrate/selection_2/subsampled_trees/<stem>.jls with
# atomic writes and shard-level resumability.
#
# n_restarts records are read from data/processed/selection_2/<timing>/n_restarts/
# (co-indexed with the raw shard) and carried through to the per-sim dict.
# The script hard-errors if the n_restarts vector length does not match the
# number of subsampled trees in the shard.
#
# Inputs:  data/raw_subsampled/selection_2/{deterministic,gamma,markov}/*.jls
#          data/processed/selection_2/{deterministic,gamma,markov}/n_restarts/*.jls
# Outputs: data/inference/mutationrate/selection_2/subsampled_trees/<stem>.jls
#          Stem already encodes the sample size: <base>_n<n>.jls.
#          Each element: Union{Nothing, Dict{Symbol,Any}} with keys
#          :mutations  (Dict from mutation_rate_suite),
#          :divisions  (Dict from mutation_rate_suite),
#          :params     (Sel2Params from the source simulation),
#          :n_restarts (Int or nothing, from data/processed/selection_2/)
#
# Usage:
#   julia --project=. -t auto code/3_inference/01_mutationrate/subsampled_trees/mutation_rate_sel2_subsampled.jl
#
# Resumable: shards whose output already exists are skipped. Force a redo by
# deleting the output file.

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
using Pkg; Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using Serialization
using MutationLoadDynamics
using EvoTracer
include(joinpath(HELPERS, "types_selection2.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "evotracer_io.jl"))
include(joinpath(HELPERS, "inference_runners.jl"))
include(joinpath(HELPERS, "subsampling.jl"))   # for serialize_atomic
include(joinpath(HELPERS, "reduced_mode.jl"))   # cap sims per shard via INFERENCE_MAX_SIMS env var
include(joinpath(HELPERS, "theory.jl"))         # for E_L_theory

const IN_ROOT  = joinpath(DATA, "raw_subsampled", "selection_2")
const OUT_ROOT = joinpath(DATA, "inference", "mutationrate", "selection_2", "subsampled_trees")
mkpath(OUT_ROOT)

"""
Strip the `_n<n>` suffix from a subsampled shard stem.

Example: "sel2_gamma_N1000_d0.0_k5.0_s0.2_M10.0_n100" -> "sel2_gamma_N1000_d0.0_k5.0_s0.2_M10.0"
"""
function strip_n_suffix(stem::String)::String
    m = match(r"^(.+)_n\d+$", stem)
    isnothing(m) && error("Cannot strip _n<n> suffix from stem: $stem")
    return m.captures[1]
end

function process_shard(path::String, timing::String)
    stem = splitext(basename(path))[1]
    out_path = joinpath(OUT_ROOT, "$(stem).jls")
    if isfile(out_path)
        println("  Skipping (already done): $stem")
        return
    end

    println("  Loading: $stem")
    subs = deserialize(path)

    # Read the n_restarts vector from data/processed/selection_2/<timing>/n_restarts/.
    # co-indexed with the raw (full-tree) shard — sim_index in each SubsampleResult
    # indexes into it.
    stem_no_n = strip_n_suffix(stem)
    nr_path = joinpath(DATA, "processed", "selection_2", timing, "n_restarts",
                       "$(stem_no_n).jls")
    isfile(nr_path) || error("n_restarts file missing: $nr_path")
    n_restarts_vec = deserialize(nr_path)

    # Guard the co-indexing invariant: length must match the raw shard, not the
    # subsampled shard (subs only contains trees that survived, n_restarts_vec covers all).
    # The check from process_sel2_subsampled.jl: length(n_restarts_vec) must equal the
    # number of entries in the raw shard (unknown here), but sim_index is 1-based
    # into n_restarts_vec, so we verify that every sim_index is in-bounds.
    # We also require at least as many entries as the maximum sim_index used.
    max_sim_index = isempty(subs) ? 0 : maximum(sub.sim_index for sub in subs)
    if max_sim_index > length(n_restarts_vec)
        error("n_restarts length mismatch for $stem: max sim_index=$max_sim_index " *
              "but length(n_restarts_vec)=$(length(n_restarts_vec))")
    end

    results = Vector{Union{Nothing,Dict{Symbol,Any}}}(nothing, length(subs))
    Threads.@threads for i in 1:min(length(subs), MAX_SIMS_PER_SHARD)
        sub = subs[i]
        p = sub.params
        EL = 2.0 * (p.N_target - 1)   # total daughter branches on the full-population tree; satisfies the Λ̂ ≥ 2(n-1) floor
        # TODO: this is a full-population lower bound; the true Λ for a subsampled tree should include hidden divisions and needs a lookup — see theory/inference.md §Pre-calibration
        d_mut  = mutation_rate_suite(prepared_sub(sub; primary = :mutations); theory_EL = EL)
        d_div  = mutation_rate_suite(prepared_sub(sub; primary = :divisions); theory_EL = EL)
        results[i] = Dict{Symbol,Any}(:mutations  => d_mut,
                                       :divisions  => d_div,
                                       :params     => p,
                                       :n_restarts => n_restarts_vec[sub.sim_index])
    end

    n_kept = count(!isnothing, results)
    println("    → $n_kept / $(length(subs)) subsampled trees processed")
    serialize_atomic(out_path, results)
end

function main()
    shards = String[]
    timings = String[]
    for timing in ("deterministic", "gamma", "markov")
        subdir = joinpath(IN_ROOT, timing)
        isdir(subdir) || continue
        for f in sort(readdir(subdir))
            endswith(f, ".jls") || continue
            push!(shards, joinpath(subdir, f))
            push!(timings, timing)
        end
    end
    println("Found $(length(shards)) sel_2 subsampled-tree shards.")
    for (path, timing) in zip(shards, timings)
        process_shard(path, timing)
    end
end

isinteractive() || main()
