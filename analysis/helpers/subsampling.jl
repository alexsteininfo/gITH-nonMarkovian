# Random subsampling of lineage trees.
#
# Requires, from the caller: `using MutationLoadDynamics` and `SubsampleResult`
# (helpers/types_subsampled.jl) plus the scenario's params types.

using AbstractTrees
using Random
using Serialization

# ── Sample sizes ──────────────────────────────────────────────────────────────

"""
    sample_sizes(N_target) -> Vector{Int}

Cell counts to draw from a tree grown to `N_target`, largest first.

An explicit table rather than `round.(fracs .* N_target)`, so the values in the
data can never drift with a change of rounding convention, and so a newly
introduced `N_target` fails loudly instead of silently acquiring sizes nobody
chose. `102 = round(0.1 * 1024)`, `1638 = round(0.1 * 16384)`,
`164 = round(0.01 * 16384)`.
"""
function sample_sizes(N_target::Int)
    N_target == 1_000  && return [100]
    N_target == 1_024  && return [102]
    N_target == 10_000 && return [1_000, 100]
    N_target == 16_384 && return [1_638, 164]
    error("no sample-size rule for N_target = $N_target — add one to sample_sizes")
end

# ── Seeding ───────────────────────────────────────────────────────────────────

"""
    sample_seed(stem, sim_index, n) -> UInt64

Seed for one draw. Depends only on values recorded in the output file, so any
single draw replays in isolation, and each threaded task builds its own rng, so
threading cannot perturb the samples.
"""
sample_seed(stem::AbstractString, sim_index::Int, n::Int) = hash((stem, sim_index, n))

# ── Induced subtree ───────────────────────────────────────────────────────────

# Copy the marked part of the tree. Left/right slots are preserved, so a node whose
# left lineage was dropped keeps `left = nothing` — the tree stays a faithful
# sub-shape of the original rather than being silently re-balanced. `NonMarkovCell`
# is immutable, so sharing `node.data` between the two trees is safe.
function _copy_marked(node::BinaryNode{NonMarkovCell},
                      marked::Set{BinaryNode{NonMarkovCell}})
    new = BinaryNode(node.data)
    if !isnothing(node.left) && node.left in marked
        new.left = _copy_marked(node.left, marked)
        new.left.parent = new
    end
    if !isnothing(node.right) && node.right in marked
        new.right = _copy_marked(node.right, marked)
        new.right.parent = new
    end
    return new
end

"""
    subsample_tree(root, n, seed) -> (new_root, sampled_ids, N_full)

Draw `n` of the tree's leaves uniformly without replacement and return the induced
lineage tree: the sampled leaves plus **every ancestor** of a sampled leaf, with the
resulting unary nodes retained rather than collapsed.

Not collapsing is what keeps the sample's observables comparable with the full
tree's. Every division ancestral to a sampled cell is still a node, so the
root-to-leaf path is unchanged and `compute_mut_per_cell` and `compute_leaf_depths`
return exactly the cell's full-tree burden and divisional depth. Collapsing would
turn `leaf_depths` into a count of bifurcations that survived sampling — a property
of the sample rather than of the cell. Because each retained edge is still one
division, `node.data.mutations` remains the per-edge mutation count and `birthtime`
differences remain the real-time edge lengths.

The founder is retained as the root even when it ends up with a single child, so
`compute_sfs` accumulates its mutations into `sfs[n]` exactly as it accumulates them
into `sfs[N]` for a full tree.

Non-destructive: `root` is left untouched, because the same tree is drawn from again
at the other sample size and the two draws are independent.
"""
function subsample_tree(root::BinaryNode{NonMarkovCell}, n::Int, seed::UInt64)
    leaves = collect(Leaves(root))
    N_full = length(leaves)
    1 <= n <= N_full ||
        error("cannot draw n = $n cells from a tree with $N_full leaves")

    # `Leaves` visits a fixed tree in a fixed order (`children` returns
    # `(left, right)`), so the draw is a pure function of (tree, seed, n).
    rng = MersenneTwister(seed)
    idx = randperm(rng, N_full)[1:n]

    # Mark each sampled leaf and its ancestors, stopping at the first node already
    # marked: total cost is the number of retained nodes, not the size of the tree.
    # `BinaryNode` is mutable, so `Set` compares by identity.
    marked = Set{BinaryNode{NonMarkovCell}}()
    for i in idx
        node = leaves[i]
        while !isnothing(node) && !(node in marked)
            push!(marked, node)
            node = node.parent
        end
    end

    new_root    = _copy_marked(root, marked)
    sampled_ids = Int64[leaves[i].data.id for i in idx]
    return new_root, sampled_ids, N_full
end
