# Equivalence gate for the migration to MutationLoadDynamics.jl v0.3.0.
#
# Originally run before deleting the local tree-statistics helpers (Task 11) and
# before removing the sampling internals of analysis/helpers/subsampling.jl
# (Task 10). Both migrations are complete; the remaining testsets continue to
# validate the package against real data:
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
# Run: julia --project=. -t auto code/data_generation/3_processing/checks/check_package_equivalence.jl

using Pkg
const ROOT = dirname(dirname(dirname(dirname(@__DIR__))))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using AbstractTrees
using Serialization
using Random
using Test

include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "types_selection2.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "subsampling.jl"))

# One small neutral shard is enough for the statistics comparison: the functions
# are tree-shape-agnostic, and N = 1000 keeps the gate to seconds rather than
# minutes. Layout is data/<stage>/<scenario>/<model>/<stem>.jls.
const STEM      = "neutral_gamma_N1000_d0.5_k5.0"
const N_SAMPLE  = 100                      # sample_sizes(1000) == [100]
const RAW_SHARD = joinpath(DATA, "raw", "neutral", "gamma",
                           "$(STEM).jls")
const SUB_SHARD = joinpath(DATA, "raw_subsampled", "neutral", "gamma",
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

# The "HR-10: tree statistics are output-preserving" testset compared the package
# functions against analysis/helpers/tree_analysis.jl. It passed on every tree of
# neutral_gamma_N1000_d0.5_k5.0.jls before that file was deleted in the v0.3.0
# migration; the comparison is no longer expressible. The remaining testsets
# compare against the arrays actually on disk, which is the stronger check anyway.

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

@testset "I-1 / HR-10: leaf_fitness matches the stored stage-2 arrays (selection 2)" begin
    # The testset above never exercises `leaf_fitness`, and its shard
    # (neutral_gamma_N1000_d0.5_k5.0) has every cell's fitness exactly 1.0 — so a
    # reordering bug in `leaf_fitness` would be invisible there even in principle.
    # `leaf_fitness` has three production call sites
    # (analysis/processing/process_sel2.jl:86,
    # analysis/processing_subsampled/process_sel1_subsampled.jl:104,
    # analysis/processing_subsampled/process_sel2_subsampled.jl:97) writing arrays
    # documented as co-indexed with `mut_per_cell`. This selection-2 shard has 200
    # trees with genuinely varying fitness, so order is observable here.
    sel2_stem = "sel2_gamma_N1000_d0.5_k5.0_s0.2_M10.0"
    sel2_raw  = joinpath(ROOT, "data", "raw", "selection_2", "gamma",
                         "$(sel2_stem).jls")
    isfile(sel2_raw) || error("missing selection-2 raw shard for the gate: $sel2_raw")

    sel2_sims  = deserialize(sel2_raw)::Vector{Sel2SimResult}
    sel2_roots = [s.tree_root for s in sel2_sims if !isnothing(s.tree_root)]
    println("  $(length(sel2_roots)) selection-2 trees")
    # Same anti-vacuity guard as the raw-shard check above: without it, a
    # corrupted or all-`nothing` shard makes the loop below a silent no-op.
    isempty(sel2_roots) &&
        error("selection-2 raw shard deserialized to zero non-nothing trees — " *
              "nothing to compare: $sel2_raw")

    # Before trusting any fitness comparison, confirm the shard actually
    # discriminates order: a shard where every tree has one fitness value would
    # pass a `leaf_fitness` check even with a reordering bug — exactly the gap
    # this testset exists to close.
    n_discriminating = count(r -> length(unique(leaf_fitness(r))) > 1, sel2_roots)
    @test n_discriminating > length(sel2_roots) ÷ 2

    checked = 0
    for quantity in ("sfs", "mut_per_cell", "leaf_depths", "leaf_fitness")
        path = joinpath(ROOT, "data", "processed", "selection_2", "gamma", quantity,
                        "$(sel2_stem).jls")
        isfile(path) || (@info "skipping (not processed): $path"; continue)
        stored = deserialize(path)
        @test length(stored) == length(sel2_roots)
        for (j, root) in enumerate(sel2_roots)
            N = length(collect(Leaves(root)))
            recomputed = quantity == "sfs"          ? sitefrequencyspectrum(root, N) :
                         quantity == "mut_per_cell" ? mutations_per_cell(root) :
                         quantity == "leaf_depths"  ? leaf_depths(root) :
                                                      leaf_fitness(root)
            @test recomputed == stored[j]
        end
        checked += 1
    end
    # Without this the testset would pass with 0 assertions if the processed
    # arrays were missing. All four quantities are confirmed present on disk in
    # this repo today.
    @test checked == 4
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
        @test mutations_per_cell(s.root) == mutations_per_cell(sub.tree_root)
        @test leaf_depths(s.root)        == leaf_depths(sub.tree_root)
        @test sitefrequencyspectrum(s.root, sub.n) == sitefrequencyspectrum(sub.tree_root, sub.n)
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
