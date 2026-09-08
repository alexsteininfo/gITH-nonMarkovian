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

# ── Exported artefacts ────────────────────────────────────────────────────────
#
# The testsets below that read these require
# `julia --project=. code/2_theory/neutral/export_trees_newick.jl` to have run.
# That split is the point: the checks are cheap, the export is not.

const NWK_DIR = joinpath(DATA, "newick", "neutral")

panel_rows() = [split(l, ',') for l in readlines(joinpath(NWK_DIR, "panels.csv"))[2:end]]

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

@testset "panels.csv describes 24 panels over 20 distinct trees" begin
    @test isfile(joinpath(NWK_DIR, "panels.csv"))
    lines = readlines(joinpath(NWK_DIR, "panels.csv"))
    @test lines[1] == "figure,model,d,sim_index,n_tips,N_pop,sampling_fraction," *
                      "t_end,root_edge,median_leaf_depth,repeated,newick_file,highlight_file"
    @test length(lines) == 25                       # header + 24 panels
    rows = panel_rows()
    @test length(unique(r[1] for r in rows)) == 4   # four figures
    @test length(unique(r[12] for r in rows)) == 20 # det d=0.9 repeats its d=0 file
    @test count(r -> r[11] == "true", rows) == 4    # one repeated panel per figure
    for r in rows
        @test isfile(joinpath(NWK_DIR, r[12]))
    end
end

@testset "every exported tree is structurally sound" begin
    for r in panel_rows()
        n_tips = parse(Int, r[5])
        s      = read(joinpath(NWK_DIR, r[12]), String)

        @test endswith(s, ";")
        # collapsed => strictly binary => n_tips - 1 internal nodes
        @test count(==('('), s) == n_tips - 1
        @test count(==(')'), s) == n_tips - 1
        @test count(==(','), s) == n_tips - 1

        labels = [m.captures[1] for m in eachmatch(r"(?:^|[(,])(\d+):", s)]
        @test length(labels) == n_tips
        @test length(unique(labels)) == n_tips
    end
end

@testset "highlight ids are a sample of the full tree's tips" begin
    for r in panel_rows()
        r[1] == "full_N1000" || continue
        @test !isempty(r[13])
        hl   = readlines(joinpath(NWK_DIR, r[13]))
        s    = read(joinpath(NWK_DIR, r[12]), String)
        tips = Set(m.captures[1] for m in eachmatch(r"(?:^|[(,])(\d+):", s))
        @test length(hl) == (r[2] == "deterministic" ? 102 : 100)
        @test all(id -> id in tips, hl)
    end
end

end
