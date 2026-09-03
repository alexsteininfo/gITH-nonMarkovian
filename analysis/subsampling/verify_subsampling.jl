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

end
