# Selection inference on full sel_1 trees.
#
# For each shard in data/raw/selection_1/<timing>/*.jls, deserialize the trees,
# convert to EvoTree (primary = :mutations), run the full selection statistic
# suite, and serialize a Vector{Union{Nothing,Dict}} to
# data/inference/selection/selection_1/full_trees/<stem>.jls with atomic writes
# and shard-level resumability.
#
# Inputs:  data/raw/selection_1/{deterministic,gamma,markov}/*.jls
# Outputs: data/inference/selection/selection_1/full_trees/<stem>.jls
#          Each element: Union{Nothing, Dict{Symbol,Any}} with keys
#          :lineages_through_time, :sibling_contrasts, :direct_contrasts,
#          :deflated_mass, :split_imbalance, :cumulative_deflation,
#          :excess_balance, :splitting_null_pvalues, :clade_ratio_score,
#          :branching_kernel_score, :empirical_contrast_null,
#          :sibling_feature_test, :lineage_fitness_scores,
#          :subsample_rank_stability, :neutral_log_deflated_mass_by_n,
#          :params     (Sel1Params from the source simulation),
#          :injection  (Sel1Injection from the source simulation)
#
# Usage:
#   julia --project=. -t auto code/3_inference/02_selection/full_trees/selection_sel1.jl
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
include(joinpath(HELPERS, "types_subsampled.jl"))  # required: evotracer_io.jl annotates SubsampleResult
include(joinpath(HELPERS, "evotracer_io.jl"))
include(joinpath(HELPERS, "inference_runners.jl"))
include(joinpath(HELPERS, "subsampling.jl"))   # for serialize_atomic
include(joinpath(HELPERS, "reduced_mode.jl"))   # cap sims per shard via INFERENCE_MAX_SIMS env var

const IN_ROOT  = joinpath(DATA, "raw", "selection_1")
const OUT_ROOT = joinpath(DATA, "inference", "selection", "selection_1", "full_trees")
mkpath(OUT_ROOT)

function process_shard(path::String)
    stem = splitext(basename(path))[1]
    out_path = joinpath(OUT_ROOT, "$(stem).jls")
    if isfile(out_path)
        println("  Skipping (already done): $stem")
        return
    end

    println("  Loading: $stem")
    sims = deserialize(path)

    results = Vector{Union{Nothing,Dict{Symbol,Any}}}(nothing, length(sims))
    Threads.@threads for i in 1:min(length(sims), MAX_SIMS_PER_SHARD)
        sim = sims[i]
        isnothing(sim.tree_root) && continue
        t = prepared_full(sim; primary = :mutations)
        results[i] = merge(selection_suite(t),
                           Dict{Symbol,Any}(:params    => sim.params,
                                            :injection => sim.injection))
    end

    n_kept = count(!isnothing, results)
    println("    → $n_kept / $(length(sims)) sims with non-nothing tree")
    serialize_atomic(out_path, results)
end

function main()
    shards = String[]
    for timing in ("deterministic", "gamma", "markov")
        subdir = joinpath(IN_ROOT, timing)
        isdir(subdir) || continue
        for f in sort(readdir(subdir))
            endswith(f, ".jls") || continue
            push!(shards, joinpath(subdir, f))
        end
    end
    println("Found $(length(shards)) sel_1 full-tree shards.")
    for path in shards
        process_shard(path)
    end
end

isinteractive() || main()
