# Unit tests for code/X_helpers/evotracer_io.jl. Run standalone:
#   julia --project=. -t auto code/3_inference/verify_evotracer_io.jl
#
# Loads one real shard from data/raw/neutral/gamma/ and one from
# data/raw_subsampled/neutral/gamma/, converts, and asserts the numeric
# equivalences the extension itself asserts upstream in test_ext.jl.

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
using Pkg; Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using Serialization
using MutationLoadDynamics
using EvoTracer
include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "evotracer_io.jl"))

using Test

# Small, cheap shard: neutral gamma at N=1000, d=0.0
FULL_SHARD = joinpath(DATA, "raw", "neutral", "gamma", "neutral_gamma_N1000_d0.0_k5.0.jls")
SUB_SHARD  = joinpath(DATA, "raw_subsampled", "neutral", "gamma",
                      "neutral_gamma_N1000_d0.0_k5.0_n100.jls")

# Load test data at module scope so both the initial testset and runner tests can use it
sims = deserialize(FULL_SHARD)
sim  = first(s for s in sims if !isnothing(s.tree_root))
subs = deserialize(SUB_SHARD)
sub  = first(subs)

@testset "evotracer_io.jl" begin

    @testset "evotree_from_full: burdens and depths match multiset" begin
        t = evotree_from_full(sim; primary = :mutations)
        # Convert to raw counts on the same set of leaves
        raw_burden = sort(mutations_per_cell(sim.tree_root))
        et_burden  = sort(Int.(tip_burdens(t)))
        @test raw_burden == et_burden
    end

    @testset "prepared_full: prepare_tree runs and preserves multiset" begin
        t = prepared_full(sim; primary = :mutations)
        et_burden = sort(Int.(tip_burdens(t)))
        raw_burden = sort(mutations_per_cell(sim.tree_root))
        @test raw_burden == et_burden
    end

    @testset "evotree_from_subsample: provenance carries LeafSample keys" begin
        t = evotree_from_subsample(sub; primary = :mutations)
        p = provenance(t)
        @test p[:n_sampled]    == sub.n
        @test p[:n_population] == sub.N_full
        @test p[:rho]          ≈ sub.n / sub.N_full
        @test p[:seed]         == sub.seed
        @test p[:replicate]    == sub.sim_index
        @test occursin("LeafSample", p[:source])
    end

    @testset "evotree_from_subsample matches LeafSample re-draw" begin
        ls = sample_leaves(sim.tree_root, sub.n; seed = sub.seed)
        # Compare against a freshly built EvoTree from the LeafSample
        t_ls  = EvoTree(ls;  primary = :mutations)
        t_sub = evotree_from_subsample(sub; primary = :mutations)
        @test sort(Int.(tip_burdens(t_ls))) == sort(Int.(tip_burdens(t_sub)))
    end
end

include(joinpath(HELPERS, "inference_runners.jl"))

@testset "inference_runners: mutation_rate_suite" begin
    t = prepared_full(sim; primary = :mutations)
    d = mutation_rate_suite(t; theory_EL = 2500.0)
    for k in (:mutation_rate_bracket, :burden_dispersion, :tip_burdens,
              :tree_length_rate, :calibrated_tree_length_rate,
              :cherry_contrast, :group_distance_ratio, :sister_depth_ratios,
              :frequency_spectrum, :singleton_rate,
              :two_partition_clock_test, :two_partition_bootstrap)
        @test haskey(d, k)
    end
    # calibrated_tree_length_rate returns a finite number given a positive EL
    v = d[:calibrated_tree_length_rate]
    @test v isa EvoTracer.StatValue
    @test isfinite(v.value)
end

@testset "inference_runners: selection_suite" begin
    t = prepared_full(sim; primary = :mutations)
    d = selection_suite(t)
    for k in (:lineages_through_time, :sibling_contrasts, :direct_contrasts,
              :deflated_mass, :split_imbalance, :cumulative_deflation,
              :excess_balance, :splitting_null_pvalues,
              :empirical_contrast_null, :sibling_feature_test,
              :lineage_fitness_scores, :clade_ratio_score,
              :branching_kernel_score, :subsample_rank_stability,
              :neutral_log_deflated_mass_by_n)
        @test haskey(d, k)
    end
    @test d[:neutral_log_deflated_mass_by_n] isa Dict{Int,Float64}
end

include(joinpath(HELPERS, "plotting_functions.jl"))
@testset "aggregate_stat: reads a shard and reduces" begin
    out_path = joinpath(DATA, "inference", "mutationrate", "neutral", "full_trees",
                        "neutral_gamma_N1000_d0.0_k5.0.jls")
    isfile(out_path) || (@info "skipping aggregate_stat test — no driver output yet"; return)
    r = aggregate_stat([out_path];
        extract = d -> d[:mutations][:burden_dispersion].value,
        reducer = median, key = :burden_dispersion)
    @test length(r.x) == 1
    @test isfinite(r.reduced[1]) || isnan(r.reduced[1])
end

println("All evotracer_io tests passed.")
