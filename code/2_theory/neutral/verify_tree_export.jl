using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using AbstractTrees
using Serialization
using Statistics
using Test

include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
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

# Re-reads the tree a panel row was built from. Deliberately duplicates the export
# driver's path construction rather than including it: a verifier that shares the
# code under test cannot catch a naming drift.
function _source_tree(r::Vector{<:AbstractString})
    model = r[2]; d = r[3]
    N = parse(Int, r[6]); i = parse(Int, r[4]); n = parse(Int, r[5])
    stem = model == "deterministic" ? "neutral_deterministic_N$(N)" :
           model == "gamma"         ? "neutral_gamma_N$(N)_d$(d)_k5.0" :
                                      "neutral_markov_N$(N)_d$(d)"
    if r[1] == "full_N1000"
        raw = deserialize(joinpath(DATA, "raw", "neutral", model, "$(stem).jls"))
        return raw[i].tree_root
    end
    subs = deserialize(joinpath(DATA, "raw_subsampled", "neutral", model, "$(stem)_n$(n).jls"))
    k = findfirst(s -> s.sim_index == i, subs)
    k === nothing && error("sim_index $i missing from $(stem)_n$(n).jls")
    return subs[k].tree_root
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

@testset "sibling cells share a birthtime" begin
    # The writer takes the left child's birthtime as the division time. If the
    # simulator ever created daughters at different times, every branch length below
    # a bifurcation would be wrong and nothing else here would notice.
    for r in panel_rows()
        r[1] == "full_N1000" || continue
        for node in PreOrderDFS(_source_tree(r))
            if node.left !== nothing && node.right !== nothing
                @test node.left.data.birthtime ≈ node.right.data.birthtime atol = 1e-12
            end
        end
    end
end

@testset "t_end and the root anchor the time axis" begin
    # Ultrametricity itself is asserted on the written file by plot_trees.R, via
    # ape::is.ultrametric — a genuinely independent parse. What can only be checked
    # here, against the source tree, is that the two ends of the axis are right.
    for r in panel_rows()
        t_end = parse(Float64, r[8])
        root  = _source_tree(r)
        # Branch lengths are emitted relative to the root's birthtime and R divides
        # by t_end, so a root born after 0 would silently shorten every panel.
        @test root.data.birthtime == 0.0
        # No cell can be born after the run stopped.
        @test tree_tmax(root) <= t_end + 1e-9
        if r[1] in ("n1000_N10000", "n100_N10000")
            # these two take t_end straight from their own sample
            @test tree_tmax(root) ≈ t_end atol = 1e-6
        end
    end
end

@testset "sampled t_end is within 1% of the population's" begin
    # The two N = 10000 figures take t_end from the sample rather than reading 1.2 GB
    # of raw trees. The gap is of order 1/(n·ln N) relative. Checked here against the
    # longest-running shard, which costs one 395 MB read — the slowest test in the file.
    raw = deserialize(joinpath(DATA, "raw", "neutral", "markov",
                               "neutral_markov_N10000_d0.9.jls"))::Vector{GrowthSimResult}
    for r in panel_rows()
        (r[2] == "markov" && r[3] == "0.9" &&
         r[1] in ("n1000_N10000", "n100_N10000")) || continue
        i      = parse(Int, r[4])
        sample = parse(Float64, r[8])
        truth  = max(tree_tmax(raw[i].tree_root), raw[i].trajectory[end].t)
        @test abs(truth - sample) / truth < 0.01
    end
end

end
