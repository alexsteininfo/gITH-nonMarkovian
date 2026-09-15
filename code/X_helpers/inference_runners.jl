# Statistic surface runners for the mutation-rate and selection inference stages.
# Kept out of the driver scripts so that adding a statistic touches one file, not
# twelve, and so the equivalence check in code/3_inference/check_inference_equivalence.jl
# has a single object to re-invoke.

using EvoTracer

# Statistics called with no kwargs. `sister_depth_ratios` is NOT here — its
# default `nperm = 999` scales as O(n_perm × n_tips) per internal node, which
# blows up on the N16384 shards (~10^11 shuffle ops for 200 sims × 2 primaries).
# Called separately below with `nperm = 0`, which returns the ratios without
# the permutation p-value — the figures only consume `.ratio`.
const MUT_STATS_NO_ARG = (:mutation_rate_bracket, :burden_dispersion, :tip_burdens,
                          :tree_length_rate, :cherry_contrast,
                          :frequency_spectrum,
                          :singleton_rate, :two_partition_clock_test,
                          :two_partition_bootstrap)

const SEL_STATS_NO_ARG = (:lineages_through_time, :sibling_contrasts,
                          :direct_contrasts, :deflated_mass, :split_imbalance,
                          :cumulative_deflation, :excess_balance,
                          :splitting_null_pvalues, :clade_ratio_score,
                          :branching_kernel_score)

_safe(f, args...; kwargs...) = try
    f(args...; kwargs...)
catch e
    e
end

"""
    mutation_rate_suite(t; theory_EL = nothing) -> Dict{Symbol,Any}

Every statistic wraps in a `try`/`catch` so one bad tree does not kill a run.
"""
function mutation_rate_suite(t::EvoTree; theory_EL::Union{Nothing,Float64} = nothing)
    out = Dict{Symbol,Any}()
    for s in MUT_STATS_NO_ARG
        out[s] = _safe(getfield(EvoTracer, s), t)
    end
    # sister_depth_ratios with nperm=0 skips the permutation null but returns the
    # per-node depth ratios that the figures actually plot. nperm=999 (default)
    # takes hours per 16k-tip tree; nperm=0 is O(n_tips).
    out[:sister_depth_ratios] = _safe(
        (tt) -> EvoTracer.sister_depth_ratios(tt; nperm = 0), t)

    # group_distance_ratio needs a region; on a lineage tree with no `region`
    # side table, pass the full tip label set (equivalent to "all tips are one
    # group") as a placeholder. The statistic degenerates to 1 there and its
    # value acts as a "did the call at least run" smoke check.
    tip_labels = EvoTracer.tiplabels(t)
    out[:group_distance_ratio] = _safe(EvoTracer.group_distance_ratio, t, tip_labels)
    if theory_EL !== nothing
        out[:calibrated_tree_length_rate] = _safe(
            EvoTracer.calibrated_tree_length_rate, t, theory_EL)
    else
        out[:calibrated_tree_length_rate] = nothing
    end
    return out
end

"""
    selection_suite(t) -> Dict{Symbol,Any}

Every statistic wraps in a `try`/`catch`. `neutral_log_deflated_mass(n)` is
precomputed once per unique clade size seen on the tree so figure code can look
them up cheaply.
"""
function selection_suite(t::EvoTree)
    out = Dict{Symbol,Any}()
    for s in SEL_STATS_NO_ARG
        out[s] = _safe(getfield(EvoTracer, s), t)
    end
    # empirical_contrast_null takes the sibling-contrast result
    if out[:sibling_contrasts] isa EvoTracer.StatTable
        out[:empirical_contrast_null] = _safe(
            EvoTracer.empirical_contrast_null, out[:sibling_contrasts])
    else
        out[:empirical_contrast_null] = nothing
    end
    # sibling_feature_test on the direct-contrast column
    out[:sibling_feature_test] = _safe(
        (t) -> EvoTracer.sibling_feature_test(t; feature = :contrast_direct), t)
    # lineage_fitness_scores: first arg must be EvoTree; pass sibling_contrasts
    # as the second arg (StatTable), with optional per-row variance from the null.
    if out[:sibling_contrasts] isa EvoTracer.StatTable
        variance_col = out[:empirical_contrast_null] isa EvoTracer.StatTable ?
            Vector{Float64}(out[:empirical_contrast_null][:variance]) : nothing
        out[:lineage_fitness_scores] = _safe(
            (tt) -> EvoTracer.lineage_fitness_scores(tt, out[:sibling_contrasts];
                                                     variance = variance_col), t)
    else
        out[:lineage_fitness_scores] = nothing
    end
    # subsample_rank_stability is only meaningful on full trees; if rho < 1 we
    # skip it. Driver decides whether to call this suite on a subsampled tree.
    out[:subsample_rank_stability] = _safe(EvoTracer.subsample_rank_stability, t)

    # neutral_log_deflated_mass by unique clade size on this tree
    sizes = if out[:sibling_contrasts] isa EvoTracer.StatTable
        Set(Int.(out[:sibling_contrasts].data.clade_size))
    else
        Set{Int}()
    end
    by_n = Dict{Int,Float64}()
    for n in sizes
        result = _safe(EvoTracer.neutral_log_deflated_mass, n)
        # neutral_log_deflated_mass returns a vector from 0 to n; store the value at the max clade size
        by_n[n] = if result isa Vector{Float64} && !isempty(result)
            result[end]
        else
            result  # store the exception if one occurred
        end
    end
    out[:neutral_log_deflated_mass_by_n] = by_n
    return out
end
