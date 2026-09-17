# Selection inference on subsampled sel_1 trees.
#
# For each shard in data/raw_subsampled/selection_1/<timing>/*.jls, deserialize
# the subsampled trees, convert to EvoTree (primary = :mutations), run the full
# selection statistic suite, and serialize a Vector{Union{Nothing,Dict}} to
# data/inference/selection/selection_1/subsampled_trees/<stem>.jls with atomic
# writes and shard-level resumability.
#
# Injection records are read from data/processed/selection_1/<timing>/injection/
# (co-indexed with the raw shard) and carried through to the per-sim dict.
# The script hard-errors if the injection vector length does not match the
# number of subsampled trees in the shard.
#
# Inputs:  data/raw_subsampled/selection_1/{deterministic,gamma,markov}/*.jls
#          data/processed/selection_1/{deterministic,gamma,markov}/injection/*.jls
# Outputs: data/inference/selection/selection_1/subsampled_trees/<stem>.jls
#          Stem already encodes the sample size: <base>_n<n>.jls.
#          Each element: Union{Nothing, Dict{Symbol,Any}} with keys
#          :lineages_through_time, :sibling_contrasts, :direct_contrasts,
#          :deflated_mass, :split_imbalance, :cumulative_deflation,
#          :excess_balance, :splitting_null_pvalues, :clade_ratio_score,
#          :branching_kernel_score, :empirical_contrast_null,
#          :sibling_feature_test, :lineage_fitness_scores,
#          :subsample_rank_stability, :neutral_log_deflated_mass_by_n,
#          :params     (Sel1Params from the source simulation),
#          :injection  (Sel1Injection or nothing, from data/processed/selection_1/)
#
# Usage:
#   julia --project=. -t auto code/3_inference/02_selection/subsampled_trees/selection_sel1_subsampled.jl
#
# Resumable: shards whose output already exists are skipped. Force a redo by
# deleting the output file.

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
using Pkg; Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using Serialization
using MutationLoadDynamics
using EvoTracer
include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "evotracer_io.jl"))
include(joinpath(HELPERS, "inference_runners.jl"))
include(joinpath(HELPERS, "subsampling.jl"))   # for serialize_atomic
include(joinpath(HELPERS, "reduced_mode.jl"))   # cap sims per shard via INFERENCE_MAX_SIMS env var

const IN_ROOT  = joinpath(DATA, "raw_subsampled", "selection_1")
const OUT_ROOT = joinpath(DATA, "inference", "selection", "selection_1", "subsampled_trees")
mkpath(OUT_ROOT)

"""
Strip the `_n<n>` suffix from a subsampled shard stem.

Example: "sel1_gamma_N1000_d0.0_k5.0_s0.4_n100" -> "sel1_gamma_N1000_d0.0_k5.0_s0.4"
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

    # Read the injection vector from data/processed/selection_1/<timing>/injection/.
    # co-indexed with the raw (full-tree) shard — sim_index in each SubsampleResult
    # indexes into it.
    stem_no_n = strip_n_suffix(stem)
    inj_path = joinpath(DATA, "processed", "selection_1", timing, "injection",
                        "$(stem_no_n).jls")
    isfile(inj_path) || error("Injection file missing: $inj_path")
    injections = deserialize(inj_path)

    # Guard the co-indexing invariant: length must match the raw shard, not the
    # subsampled shard (subs only contains trees that survived, injections covers all).
    # Verify that every sim_index is in-bounds.
    max_sim_index = isempty(subs) ? 0 : maximum(sub.sim_index for sub in subs)
    if max_sim_index > length(injections)
        error("Injection length mismatch for $stem: max sim_index=$max_sim_index " *
              "but length(injections)=$(length(injections))")
    end

    results = Vector{Union{Nothing,Dict{Symbol,Any}}}(nothing, length(subs))
    Threads.@threads for i in 1:min(length(subs), MAX_SIMS_PER_SHARD)
        sub = subs[i]
        t = prepared_sub(sub; primary = :mutations)
        results[i] = merge(selection_suite(t),
                           Dict{Symbol,Any}(:params    => sub.params,
                                            :injection => injections[sub.sim_index]))
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
            shard_should_skip(f, true) && continue
            push!(shards, joinpath(subdir, f))
            push!(timings, timing)
        end
    end
    println("Found $(length(shards)) sel_1 subsampled-tree shards.")
    for (path, timing) in zip(shards, timings)
        process_shard(path, timing)
    end
end

isinteractive() || main()
