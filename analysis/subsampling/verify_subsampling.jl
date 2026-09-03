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

@testset "serialize_atomic writes atomically and leaves no tmp file" begin
    dir  = mktempdir()
    path = joinpath(dir, "demo.jls")
    tmp  = path * ".tmp"

    @test !isfile(path)
    serialize_atomic(path, [1, 2, 3])
    @test isfile(path)
    @test !isfile(tmp)
    @test deserialize(path) == [1, 2, 3]

    # A forced redo — the pipeline's resumability idiom is "delete the output, rerun"
    # — still leaves the destination complete and no tmp file behind.
    serialize_atomic(path, [9, 9])
    @test isfile(path)
    @test !isfile(tmp)
    @test deserialize(path) == [9, 9]

    rm(dir; recursive = true)
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

        # Same seed reproduces a draw (above); different seeds must actually give a
        # different draw. Seeds 42 and 2 were confirmed by hand to land on different
        # pairs of the fixture's 3 leaves ([3, 5] vs [3, 4]), so this cannot flake on
        # an unlucky coincidence of seed and fixture.
        c, ids_c, _ = subsample_tree(full, 2, UInt64(2))
        @test sort(ids_a) == [3, 5]
        @test sort(ids_c) == [3, 4]
        @test sort(ids_a) != sort(ids_c)

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

# ── Helpers: id-labelled reference values from a full tree ────────────────────
# `compute_leaf_depths` returns an unlabelled vector in its own traversal order, so
# comparing a subsample cell-by-cell needs id keys. `id_burden_map` reproduces
# `mutations_per_cell`: it sums every ancestor's mutations including the root's.

function id_depth_map(root::BinaryNode{NonMarkovCell})
    m = Dict{Int64, Int}()
    stack = Tuple{BinaryNode{NonMarkovCell}, Int}[(root, 0)]
    while !isempty(stack)
        node, d = pop!(stack)
        if isnothing(node.left) && isnothing(node.right)
            m[node.data.id] = d
        else
            isnothing(node.left)  || push!(stack, (node.left,  d + 1))
            isnothing(node.right) || push!(stack, (node.right, d + 1))
        end
    end
    return m
end

function id_burden_map(root::BinaryNode{NonMarkovCell})
    m = Dict{Int64, Int}()
    for leaf in Leaves(root)
        muts = leaf.data.mutations
        node = leaf
        while !isnothing(node.parent)
            node = node.parent
            muts += node.data.mutations
        end
        m[leaf.data.id] = muts
    end
    return m
end

@testset "real shard against data/processed" begin
    stem     = "neutral_gamma_N1000_d0.5_k5.0"
    raw_path = joinpath(ROOT, "data", "raw", "neutral", "gamma", stem * ".jls")
    proc_dir = joinpath(ROOT, "data", "processed", "neutral", "gamma")

    if !isfile(raw_path) || !isdir(proc_dir)
        @info "skipping real-shard tests: $raw_path or $proc_dir not present"
    else
        sims       = deserialize(raw_path)::Vector{GrowthSimResult}
        ref_sfs    = deserialize(joinpath(proc_dir, "sfs",          stem * ".jls"))
        ref_mpc    = deserialize(joinpath(proc_dir, "mut_per_cell", stem * ".jls"))
        ref_depths = deserialize(joinpath(proc_dir, "leaf_depths",  stem * ".jls"))

        # `process_neutral.jl` drops nothing-tree sims in source order; do the same
        # so index i here is index i there.
        kept = [i for i in eachindex(sims) if !isnothing(sims[i].tree_root)]
        @test length(kept) == length(ref_sfs)

        @testset "n = N_full reproduces the processed arrays exactly" begin
            for j in 1:3
                i    = kept[j]
                root = sims[i].tree_root
                sub, ids, N_full = subsample_tree(root, 1_000, UInt64(j))
                @test N_full == 1_000
                @test sort(ids) == sort([l.data.id for l in Leaves(root)])
                @test compute_sfs(sub, N_full)  == ref_sfs[j]
                @test compute_mut_per_cell(sub) == ref_mpc[j]
                @test compute_leaf_depths(sub)  == ref_depths[j]
            end
        end

        @testset "n = 100 leaves cell-level quantities unchanged" begin
            for j in 1:3
                i     = kept[j]
                root  = sims[i].tree_root
                depth = id_depth_map(root)
                burden = id_burden_map(root)

                sub, ids, N_full = subsample_tree(root, 100, sample_seed(stem, i, 100))
                @test N_full == 1_000
                @test length(ids) == 100
                @test length(unique(ids)) == 100                  # without replacement

                sub_ids = Int64[l.data.id for l in Leaves(sub)]
                @test sort(sub_ids) == sort(ids)                  # leaves are the draw
                @test all(haskey(depth, id) for id in sub_ids)    # drawn from this tree

                # id-matched equality with the full tree, both quantities
                @test id_depth_map(sub)  == Dict(id => depth[id]  for id in sub_ids)
                @test id_burden_map(sub) == Dict(id => burden[id] for id in sub_ids)
                @test compute_mut_per_cell(sub) == Int[burden[id] for id in sub_ids]

                # SFS bookkeeping
                sfs = compute_sfs(sub, 100)
                @test length(sfs) == 100
                @test sum(k * sfs[k] for k in 1:100) == sum(compute_mut_per_cell(sub))
                @test sum(sfs) == sum(nd.data.mutations for nd in PreOrderDFS(sub))
                @test sum(sfs) < sum(ref_sfs[j])                  # sampling loses mutations

                # structure
                @test isnothing(sub.parent)
                @test sub.data.id == root.data.id                 # founder retained

                # the source tree is untouched
                @test compute_sfs(root, 1_000) == ref_sfs[j]
            end
        end

        @testset "subsample_shard preserves order and filters nothing trees" begin
            res = subsample_shard(sims[1:5], stem, 100)
            @test res isa Vector{SubsampleResult{SimParams}}
            @test length(res) == count(!isnothing, [s.tree_root for s in sims[1:5]])
            @test [r.sim_index for r in res] == [i for i in 1:5 if !isnothing(sims[i].tree_root)]
            for r in res
                @test r.n      == 100
                @test r.N_full == 1_000
                @test r.seed   === sample_seed(stem, r.sim_index, 100)
                @test r.params == sims[r.sim_index].params
                @test length(r.sampled_ids) == 100
                @test sort([l.data.id for l in Leaves(r.tree_root)]) == sort(r.sampled_ids)
            end
        end

        @testset "subsample_file writes one size and is resumable" begin
            out = mktempdir()
            small = joinpath(out, "src")
            mkpath(small)
            serialize(joinpath(small, stem * ".jls"), sims[1:4])

            subsample_file(joinpath(small, stem * ".jls"), out, [100])
            f = joinpath(out, stem * "_n100.jls")
            @test isfile(f)
            first_ids = [r.sampled_ids for r in deserialize(f)]
            @test length(first_ids) == 4

            mtime_before = mtime(f)
            subsample_file(joinpath(small, stem * ".jls"), out, [100])   # resumes: no work
            @test mtime(f) == mtime_before

            rm(f)
            subsample_file(joinpath(small, stem * ".jls"), out, [100])   # redo after delete
            @test [r.sampled_ids for r in deserialize(f)] == first_ids   # same draws
            rm(out; recursive = true)
        end

        @testset "subsample_file skips only the sizes already present, regenerating the rest" begin
            out = mktempdir()
            small = joinpath(out, "src")
            mkpath(small)
            serialize(joinpath(small, stem * ".jls"), sims[1:4])

            subsample_file(joinpath(small, stem * ".jls"), out, [100, 50])
            f100 = joinpath(out, stem * "_n100.jls")
            f50  = joinpath(out, stem * "_n50.jls")
            @test isfile(f100)
            @test isfile(f50)
            ids100_before = [r.sampled_ids for r in deserialize(f100)]
            ids50_before  = [r.sampled_ids for r in deserialize(f50)]

            rm(f50)
            mtime100_before = mtime(f100)
            subsample_file(joinpath(small, stem * ".jls"), out, [100, 50])   # only n=50 redone

            @test mtime(f100) == mtime100_before                              # survivor untouched
            @test [r.sampled_ids for r in deserialize(f100)] == ids100_before # survivor's content unchanged
            @test isfile(f50)
            @test [r.sampled_ids for r in deserialize(f50)] == ids50_before   # regenerated with same draws
            rm(out; recursive = true)
        end
    end
end

@testset "subsample_shard filters synthetic nothing-tree holes, preserving order and co-indexing" begin
    # The real shard used above is stipulated to have zero nothing-tree entries, so
    # it can never exercise the nothing-skip / co-indexing branch of subsample_shard.
    # Build a synthetic shard with real holes instead: nothing at the first, an
    # interior, and the last position, with kept sims sandwiched between them and
    # each carrying distinguishable params so the sim_index → params mapping is
    # actually checked, not just the output length.
    p(nu) = SimParams(1.0, 0.5, 5.0, nu, 0.0, 1_000, NaN)
    sims_holes = GrowthSimResult[
        GrowthSimResult(TrajectoryPoint[], nothing,       p(1.0)),  # 1: hole (first)
        GrowthSimResult(TrajectoryPoint[], fixture_tree(), p(2.0)), # 2: kept
        GrowthSimResult(TrajectoryPoint[], nothing,       p(3.0)),  # 3: hole (interior)
        GrowthSimResult(TrajectoryPoint[], fixture_tree(), p(4.0)), # 4: kept
        GrowthSimResult(TrajectoryPoint[], nothing,       p(5.0)),  # 5: hole (last)
    ]

    res = subsample_shard(sims_holes, "holes-stem", 2)
    @test res isa Vector{SubsampleResult{SimParams}}
    @test length(res) == 2
    @test [r.sim_index for r in res] == [2, 4]          # ascending source order, holes removed
    for r in res
        @test r.params == sims_holes[r.sim_index].params   # mapped back to its own source sim
    end
end

end
