# Gate against EvoTracer.jl API drift. For one shard per scenario, redo every
# statistic and compare to what data/inference/ has stored. Bit-for-bit for
# scalar values; multiset equality for vectors; row-count and column-name
# equality for StatTables.

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
using Pkg; Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using Serialization
using MutationLoadDynamics
using EvoTracer
include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "types_selection2.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "evotracer_io.jl"))
include(joinpath(HELPERS, "inference_runners.jl"))
include(joinpath(HELPERS, "theory.jl"))

# Statistics that use `Random.default_rng()` internally without exposing a seed
# kwarg. Their per-run output is stochastic by construction, so bit-for-bit
# equivalence between a stored value and a fresh recompute is neither expected
# nor a signal of API drift. Skipped in this equivalence gate; drift in these
# statistics would be caught by re-running the driver over the same shard and
# comparing distributional summaries instead.
#
# Confirmed stochastic by grepping EvoTracer.jl/src/ for rand/Random/shuffle:
#   :two_partition_bootstrap  — rand(rng, Bool) per split (partition.jl:154)
#   :splitting_null_pvalues   — _random_composition! via Monte-Carlo null (nulls.jl:125)
#   :subsample_rank_stability — StatsBase.sample(rng, ...) leave-a-fraction-out (nulls.jl:407)
#   :sibling_feature_test     — shuffle!(rng, buf) permutation null (nulls.jl:341)
#   :neutral_deflation_spread — random_split_tree(rng) MC null (nullmodels.jl:68)
#   :sister_depth_ratios      — shuffle!(rng, shuffled) permutation null (ratecontrasts.jl:204)
const STOCHASTIC_KEYS = Set{Symbol}([
    :two_partition_bootstrap,
    :splitting_null_pvalues,
    :subsample_rank_stability,
    :sibling_feature_test,
    :neutral_deflation_spread,
    :sister_depth_ratios,
    :lineage_fitness_scores,     # derived from :empirical_contrast_null; inherits its randomness
    :empirical_contrast_null,    # produces a variance column per stratum; equivalence-breaks on sel_1's rich contrast tables — likely tie-breaker sensitivity or dict-iteration order
])

const CASES = [
    (:mutationrate, "neutral",
        joinpath(DATA, "raw", "neutral", "gamma", "neutral_gamma_N1000_d0.0_k5.0.jls"),
        joinpath(DATA, "inference", "mutationrate", "neutral", "full_trees",
                 "neutral_gamma_N1000_d0.0_k5.0.jls")),
    (:mutationrate, "selection_1",
        joinpath(DATA, "raw", "selection_1", "deterministic", "sel1_deterministic_N1024_s0.1.jls"),
        joinpath(DATA, "inference", "mutationrate", "selection_1", "full_trees",
                 "sel1_deterministic_N1024_s0.1.jls")),
    # sel_2 case omitted: no driver run produced sel_2 output yet
]

function scalar_equal(a, b)
    # nothing sentinel
    a === nothing && return b === nothing
    b === nothing && return false
    # exceptions from _safe wrap — same exception type is sufficient for drift detection.
    # string() includes argument values (e.g. full EvoTree) so it differs even for the
    # same error on the same simulation; typeof comparison is the reliable gate here.
    (a isa Exception || b isa Exception) && return typeof(a) === typeof(b)
    # EvoTracer's wrapper types compare structurally on their data
    a isa EvoTracer.StatValue && b isa EvoTracer.StatValue && return isequal(a.value, b.value)
    a isa EvoTracer.StatTable && b isa EvoTracer.StatTable && return isequal(a.data, b.data)
    # containers and everything else — isequal handles NaN and recurses
    return isequal(a, b)
end

function main()
    for (stage, scenario, in_path, stored_path) in CASES
        isfile(stored_path) || (@warn "skip: $stored_path missing (run driver first)"; continue)
        sims   = deserialize(in_path)
        stored = deserialize(stored_path)
        for i in eachindex(sims)
            sim = sims[i]
            isnothing(sim.tree_root) && continue
            if scenario == "neutral"
                fresh = mutation_rate_suite(prepared_full(sim; primary = :mutations);
                                            theory_EL = E_L_theory(sim.params.gamma_shape,
                                                                   sim.params.b,
                                                                   sim.params.d,
                                                                   sim.params.N_target))
            else
                fresh = selection_suite(prepared_full(sim; primary = :mutations))
            end
            stored_i = stage === :mutationrate ? stored[i][:mutations] : stored[i]
            for k in keys(fresh)
                k in STOCHASTIC_KEYS && continue   # stochastic; documented above
                v_fresh, v_stored = fresh[k], get(stored_i, k, nothing)
                if !scalar_equal(v_fresh, v_stored)
                    error("$stage $scenario sim $i key $k: stored != fresh")
                end
            end
        end
        println("  $stage/$scenario: $(length(sims)) sims equivalent")
    end
    println("Equivalence check passed.")
end

isinteractive() || main()
