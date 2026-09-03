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

# ── Fixture ───────────────────────────────────────────────────────────────────
# A hand-built tree whose observables are known by hand and were confirmed
# against the real helpers:
#
#   root (id 1, mut 5, t 0.0)
#   ├── L  (id 2, mut 1, t 1.0)
#   │   ├── LL (id 4, mut 2, t 2.0)   leaf
#   │   └── LR (id 5, mut 3, t 2.1)   leaf
#   └── R  (id 3, mut 7, t 1.2)       leaf
#
#   Leaves(root) order = [4, 5, 3]
#   compute_mut_per_cell = [8, 9, 12]   (root's 5 is included in every burden)
#   compute_leaf_depths  = [1, 2, 2]    (own stack order, not Leaves order)
#   compute_sfs(root, 3) = [12, 1, 5]

function fixture_tree()
    root = BinaryNode(NonMarkovCell(1, 0.0, 5, 1.0))
    L = leftchild!(root,  NonMarkovCell(2, 1.0, 1, 1.0))
    rightchild!(root,     NonMarkovCell(3, 1.2, 7, 1.0))
    leftchild!(L,         NonMarkovCell(4, 2.0, 2, 1.0))
    rightchild!(L,        NonMarkovCell(5, 2.1, 3, 1.0))
    return root
end

@testset "subsampling" begin

@testset "fixture matches the helpers" begin
    root = fixture_tree()
    @test [l.data.id for l in Leaves(root)] == [4, 5, 3]
    @test compute_mut_per_cell(root)        == [8, 9, 12]
    @test compute_leaf_depths(root)         == [1, 2, 2]
    @test compute_sfs(root, 3)              == [12, 1, 5]
end

@testset "SubsampleResult round-trips through Serialization" begin
    root = fixture_tree()

    neutral = SubsampleResult(root, SimParams(1.0, 0.5, 5.0, 2.0, 0.0, 1_000, NaN),
                              3, 3, 7, UInt64(0xabc), Int64[4, 5, 3])
    sel1    = SubsampleResult(root, Sel1Params(1.0, 0.5, 5.0, 2.0, 1_000, :gamma, 0.3, 44, 2),
                              3, 3, 7, UInt64(0xabc), Int64[4, 5, 3])
    sel2    = SubsampleResult(root, Sel2Params(1.0, 0.5, 5.0, 1.0, 1_000, :gamma, 0.2, 10.0, 1.0, 3, UInt64(9)),
                              3, 3, 7, UInt64(0xabc), Int64[4, 5, 3])

    for (r, P) in ((neutral, SimParams), (sel1, Sel1Params), (sel2, Sel2Params))
        @test r isa SubsampleResult{P}
        path = tempname() * ".jls"
        serialize(path, [r])
        back = deserialize(path)::Vector{SubsampleResult{P}}
        @test length(back) == 1
        @test back[1].params      == r.params
        @test back[1].n           == 3
        @test back[1].N_full      == 3
        @test back[1].sim_index   == 7
        @test back[1].seed        === UInt64(0xabc)
        @test back[1].sampled_ids == Int64[4, 5, 3]
        # the tree survives the round trip, parent links included
        @test compute_sfs(back[1].tree_root, 3) == [12, 1, 5]
        @test compute_mut_per_cell(back[1].tree_root) == [8, 9, 12]
        rm(path)
    end
end

@testset "sample_sizes is an exact table" begin
    @test sample_sizes(1_000)  == [100]
    @test sample_sizes(1_024)  == [102]
    @test sample_sizes(10_000) == [1_000, 100]
    @test sample_sizes(16_384) == [1_638, 164]
    @test_throws ErrorException sample_sizes(5_000)
    # largest first, and every size a strict subsample
    for N in (1_000, 1_024, 10_000, 16_384)
        ns = sample_sizes(N)
        @test issorted(ns; rev = true)
        @test all(0 .< ns .< N)
    end
end

@testset "sample_seed is deterministic and draw-specific" begin
    @test sample_seed("stem", 1, 100) === sample_seed("stem", 1, 100)
    @test sample_seed("stem", 1, 100) !== sample_seed("stem", 2, 100)
    @test sample_seed("stem", 1, 100) !== sample_seed("stem", 1, 164)
    @test sample_seed("stem", 1, 100) !== sample_seed("other", 1, 100)
    @test sample_seed("stem", 1, 100) isa UInt64
end

@testset "subsample_tree on the fixture" begin
    full = fixture_tree()
    full_sfs    = compute_sfs(full, 3)
    full_mpc    = compute_mut_per_cell(full)
    full_depths = compute_leaf_depths(full)

    # id-labelled reference values from the full tree
    depth_of  = Dict(4 => 2, 5 => 2, 3 => 1)
    burden_of = Dict(4 => 8, 5 => 9, 3 => 12)

    @testset "n = N_full reproduces the source tree exactly" begin
        sub, ids, N_full = subsample_tree(full, 3, UInt64(1))
        @test N_full == 3
        @test sort(ids) == [3, 4, 5]
        @test compute_sfs(sub, 3)       == full_sfs
        @test compute_mut_per_cell(sub) == full_mpc
        @test compute_leaf_depths(sub)  == full_depths
        @test [l.data.id for l in Leaves(sub)] == [4, 5, 3]
    end

    @testset "single-cell samples over many seeds" begin
        drawn = Int64[]
        for s in UInt64(1):UInt64(30)
            sub, ids, N_full = subsample_tree(full, 1, s)
            @test N_full == 3
            @test length(ids) == 1
            append!(drawn, ids)

            # Leaves are exactly the draw. This is also what rules out a retained
            # node with no sampled descendant: such a node would surface here as a
            # leaf whose id is not in `ids`.
            leaf_ids = [l.data.id for l in Leaves(sub)]
            @test leaf_ids == ids
            @test isnothing(sub.parent)                 # new root is detached
            @test sub.data.id == 1                      # founder retained
            # cell-level quantities are unchanged by sampling
            @test compute_leaf_depths(sub)  == [depth_of[ids[1]]]
            @test compute_mut_per_cell(sub) == [burden_of[ids[1]]]
            # SFS bookkeeping: every mutation counted once per carrier, both sides
            sfs = compute_sfs(sub, 1)
            @test length(sfs) == 1
            @test sum(k * sfs[k] for k in 1:1) == sum(compute_mut_per_cell(sub))
        end
        @test sort(unique(drawn)) == [3, 4, 5]           # all leaves reachable
    end

    @testset "two-cell samples keep the branching node" begin
        # whenever both sampled cells sit under L, sfs[2] must carry L's mutation
        seen = false
        for s in UInt64(1):UInt64(60)
            sub, ids, _ = subsample_tree(full, 2, s)
            @test sort([l.data.id for l in Leaves(sub)]) == sort(ids)
            sfs = compute_sfs(sub, 2)
            @test sum(k * sfs[k] for k in 1:2) == sum(compute_mut_per_cell(sub))
            if sort(ids) == [4, 5]
                seen = true
                @test sfs == [5, 6]   # sfs[1] = 2 + 3, sfs[2] = L's 1 + root's 5
            end
        end
        @test seen
    end

    @testset "determinism and bounds" begin
        a, ids_a, _ = subsample_tree(full, 2, UInt64(42))
        b, ids_b, _ = subsample_tree(full, 2, UInt64(42))
        @test ids_a == ids_b
        @test compute_sfs(a, 2) == compute_sfs(b, 2)
        @test_throws ErrorException subsample_tree(full, 4, UInt64(1))
        @test_throws ErrorException subsample_tree(full, 0, UInt64(1))
    end

    @testset "the source tree is not mutated" begin
        subsample_tree(full, 1, UInt64(3))
        subsample_tree(full, 2, UInt64(4))
        @test compute_sfs(full, 3)       == full_sfs
        @test compute_mut_per_cell(full) == full_mpc
        @test compute_leaf_depths(full)  == full_depths
    end

    @testset "parent links in the rebuilt tree" begin
        sub, _, _ = subsample_tree(full, 2, UInt64(7))
        @test isnothing(sub.parent)
        for node in PreOrderDFS(sub)
            isnothing(node.left)  || @test node.left.parent  === node
            isnothing(node.right) || @test node.right.parent === node
        end
    end
end

end
