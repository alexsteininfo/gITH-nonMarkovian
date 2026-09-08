using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using AbstractTrees
using Statistics
using Test

include(joinpath(HELPERS, "newick.jl"))

# ── Fixture ───────────────────────────────────────────────────────────────────
#
# A hand-built tree whose Newick form is known by hand. Unlike the fixture in
# verify_subsampling.jl this one is *physical*: both daughters of a division share
# a birthtime, which is what the simulator produces and what the writer assumes.
#
#   root id1  born 0.0, divides at 1.0
#     L  id2  born 1.0, divides at 2.5      binary
#       id4  born 2.5   leaf
#       id5  born 2.5   leaf
#     R  id3  born 1.0, divides at 1.8      unary — the other daughter died
#       id6  born 1.8   leaf
#   tmax = 3.0
#
# Collapsed, R contributes its 0.8 lifetime to id6's edge: 0.8 + (3.0 - 1.8) = 2.0.

function fixture_tree()
    root = BinaryNode(NonMarkovCell(1, 0.0, 0, 1.0))
    L = leftchild!(root,  NonMarkovCell(2, 1.0, 0, 1.0))
    R = rightchild!(root, NonMarkovCell(3, 1.0, 0, 1.0))
    leftchild!(L,  NonMarkovCell(4, 2.5, 0, 1.0))
    rightchild!(L, NonMarkovCell(5, 2.5, 0, 1.0))
    leftchild!(R,  NonMarkovCell(6, 1.8, 0, 1.0))
    return root
end

@testset "tree export" begin

@testset "newick_string on the fixture" begin
    root = fixture_tree()
    @test newick_string(root; tmax = 3.0) == "((4:0.5,5:0.5):1.5,6:2):1;"
end

@testset "tree_tmax is the latest leaf birth" begin
    @test tree_tmax(fixture_tree()) == 2.5
end

@testset "the fixture tree is ultrametric under the writer's convention" begin
    # every root-to-leaf branch-length sum must equal tmax
    s = newick_string(fixture_tree(); tmax = 3.0)
    @test s == "((4:0.5,5:0.5):1.5,6:2):1;"   # 1+1.5+0.5 == 1+2 == 3.0
end

@testset "write_atomic leaves no .tmp behind" begin
    path = tempname() * ".nwk"
    write_atomic(path, "hello")
    @test read(path, String) == "hello"
    @test !isfile(path * ".tmp")
    rm(path)
end

@testset "write_newick round-trips the string" begin
    path = tempname() * ".nwk"
    write_newick(path, fixture_tree(); tmax = 3.0)
    @test read(path, String) == "((4:0.5,5:0.5):1.5,6:2):1;"
    rm(path)
end

end
