# Equivalence gate for the migration to MutationLoadDynamics.jl v0.3.0.
#
# Run this BEFORE deleting analysis/helpers/tree_analysis.jl or the sampling
# internals of analysis/helpers/subsampling.jl. It answers two questions against
# real data:
#
#   HR-10  do the package's tree statistics return exactly what the local helpers
#          return, element for element and in the same order? 4.3 GB of stage-2
#          arrays were produced by the local versions, and stage 2 is resumable by
#          `isfile`, so a reordered output would mix conventions inside one
#          dataset with no error anywhere.
#   HR-1   does the package sampler reproduce the draws already serialized under
#          data/raw_subsampled/? 4.2 GB of sampled trees and every figure derived
#          from them depend on that.
#
# Exits non-zero on any mismatch. A mismatch means Phase A is wrong; it never
# means the data should be regenerated.
#
# Run: julia --project=. analysis/checks/check_package_equivalence.jl

using Pkg
Pkg.activate(dirname(dirname(@__DIR__)))

using MutationLoadDynamics
using AbstractTrees
using Serialization
using Random
using Test

const ROOT    = dirname(dirname(@__DIR__))
const HELPERS = joinpath(dirname(@__DIR__), "helpers")

include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "types_selection2.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "tree_analysis.jl"))
include(joinpath(HELPERS, "subsampling.jl"))

# One small neutral shard is enough for the statistics comparison: the functions
# are tree-shape-agnostic, and N = 1000 keeps the gate to seconds rather than
# minutes. Layout is data/<stage>/<scenario>/<model>/<stem>.jls.
const STEM      = "neutral_gamma_N1000_d0.5_k5.0"
const N_SAMPLE  = 100                      # sample_sizes(1000) == [100]
const RAW_SHARD = joinpath(ROOT, "data", "raw", "neutral", "gamma",
                           "$(STEM).jls")
const SUB_SHARD = joinpath(ROOT, "data", "raw_subsampled", "neutral", "gamma",
                           "$(STEM)_n$(N_SAMPLE).jls")

isfile(RAW_SHARD) || error("missing raw shard for the gate: $RAW_SHARD")
isfile(SUB_SHARD) || error("missing subsampled shard for the gate: $SUB_SHARD")

println("Deserializing $(basename(RAW_SHARD)) …")
sims = deserialize(RAW_SHARD)
roots = [s.tree_root for s in sims if !isnothing(s.tree_root)]
println("  $(length(roots)) trees")
# Without this, a corrupted or all-`nothing` raw shard makes every loop below a
# silent no-op: the testsets that iterate `roots` would report 0/0 and pass,
# and the gate would print "safe to migrate" having compared nothing.
isempty(roots) &&
    error("raw shard deserialized to zero non-nothing trees — nothing to compare: $RAW_SHARD")

@testset "package equivalence gate" begin

@testset "HR-10: tree statistics are output-preserving" begin
    for (i, root) in enumerate(roots)
        N = length(collect(Leaves(root)))
        # Tautological by construction: compute_mut_per_cell is a one-line
        # pass-through to this very function, and mutations_per_cell pre-dates
        # this migration. Kept for symmetry with the other five. The real
        # coverage for this quantity is the stored-array comparison below and
        # the induced-tree check in the HR-1 testset.
        @test mutations_per_cell(root)          == compute_mut_per_cell(root)
        @test leaf_depths(root)                 == compute_leaf_depths(root)
        @test sitefrequencyspectrum(root, N)    == compute_sfs(root, N)
        @test branch_spectrum(root, N)          == compute_branch_spectrum(root, N)
        @test leaf_fitness(root)                == compute_leaf_fitness(root)
        @test filtered_mutations_per_cell(root, 0.3) ==
              compute_filtered_mut_per_cell(root, 0.3)
    end
end

@testset "HR-10: package statistics match the stored stage-2 arrays" begin
    # The strongest form of the check: compare against what is actually on disk,
    # not just against the helper that produced it.
    checked = 0
    for quantity in ("sfs", "mut_per_cell", "leaf_depths")
        path = joinpath(ROOT, "data", "processed", "neutral", "gamma", quantity,
                        "$(STEM).jls")
        isfile(path) || (@info "skipping (not processed): $path"; continue)
        stored = deserialize(path)
        @test length(stored) == length(roots)
        for (j, root) in enumerate(roots)
            N = length(collect(Leaves(root)))
            recomputed = quantity == "sfs"          ? sitefrequencyspectrum(root, N) :
                         quantity == "mut_per_cell" ? mutations_per_cell(root) :
                                                      leaf_depths(root)
            @test recomputed == stored[j]
        end
        checked += 1
    end
    # Without this the testset passes with 0 assertions when the processed arrays
    # are missing, and the gate prints "safe to migrate" having verified nothing.
    # All three quantities are confirmed present on disk in this repo today.
    @test checked == 3
end

@testset "HR-1: the package sampler reproduces the stored draws" begin
    stored = deserialize(SUB_SHARD)
    @test !isempty(stored)
    for sub in stored
        root = sims[sub.sim_index].tree_root
        @test !isnothing(root)
        s = sample_leaves(root, sub.n; seed = sub.seed)
        @test s.sampled_ids == sub.sampled_ids
        @test s.N_full      == sub.N_full
        # And the induced trees agree on every observable.
        @test mutations_per_cell(s.root) == compute_mut_per_cell(sub.tree_root)
        @test leaf_depths(s.root)        == compute_leaf_depths(sub.tree_root)
        @test sitefrequencyspectrum(s.root, sub.n) == compute_sfs(sub.tree_root, sub.n)
    end
end

@testset "HR-2: the study seed derivation still reaches the same draws" begin
    stored = deserialize(SUB_SHARD)
    for sub in stored
        # subsample_file derives the seed from the *raw shard's* stem.
        @test sub.seed == sample_seed(STEM, sub.sim_index, sub.n)
        root = sims[sub.sim_index].tree_root
        s = sample_leaves(root, sub.n; seed = sample_seed(STEM, sub.sim_index, sub.n))
        @test s.sampled_ids == sub.sampled_ids
    end
end

end

println("\nGate passed — safe to migrate.")
