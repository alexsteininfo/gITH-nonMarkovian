# Tree-based analysis helpers for MutationLoadDynamics.jl simulations.
# Requires: using MutationLoadDynamics, using AbstractTrees

using AbstractTrees

# ── Mutation counts ────────────────────────────────────────────────────────────

function compute_mut_per_cell(root::BinaryNode{NonMarkovCell})
    return mutations_per_cell(root)
end

# ── Site frequency spectrum ────────────────────────────────────────────────────

function _fill_sfs!(node::BinaryNode{NonMarkovCell}, sfs::Vector{Int64})
    if isnothing(node.left) && isnothing(node.right)
        sfs[1] += node.data.mutations
        return 1
    end
    count = 0
    isnothing(node.left)  || (count += _fill_sfs!(node.left,  sfs))
    isnothing(node.right) || (count += _fill_sfs!(node.right, sfs))
    count > 0 && (sfs[count] += node.data.mutations)
    return count
end

function compute_sfs(root::BinaryNode{NonMarkovCell}, N::Int)
    sfs = zeros(Int64, N)
    _fill_sfs!(root, sfs)
    return sfs
end

# ── Leaf divisional depths ─────────────────────────────────────────────────────
# Returns the number of division events from root to each leaf.
# For neutral simulations (ν=0), this is the primary quantity: mutations per cell
# are obtained post-hoc by drawing Poisson(m * depth) for each leaf.

function compute_leaf_depths(root::BinaryNode{NonMarkovCell})
    depths = Int[]
    stack  = Tuple{BinaryNode{NonMarkovCell}, Int}[(root, 0)]
    while !isempty(stack)
        node, d = pop!(stack)
        if isnothing(node.left) && isnothing(node.right)
            push!(depths, d)
        else
            isnothing(node.left)  || push!(stack, (node.left,  d + 1))
            isnothing(node.right) || push!(stack, (node.right, d + 1))
        end
    end
    return depths
end

# ── Topological site frequency spectrum ────────────────────────────────────────
# branch_spectrum[n] = number of internal nodes subtending exactly n leaves.
# For neutral mutations (ν=0), this encodes the full tree topology: the expected
# SFS under any neutral mutation model m is E[ξ_n] = m * branch_spectrum[n].

function _fill_branch_spectrum!(node::BinaryNode{NonMarkovCell}, bs::Vector{Int})
    isnothing(node.left) && isnothing(node.right) && return 1
    count = 0
    isnothing(node.left)  || (count += _fill_branch_spectrum!(node.left,  bs))
    isnothing(node.right) || (count += _fill_branch_spectrum!(node.right, bs))
    count > 0 && (bs[count] += 1)
    return count
end

function compute_branch_spectrum(root::BinaryNode{NonMarkovCell}, N::Int)
    bs = zeros(Int, N)
    _fill_branch_spectrum!(root, bs)
    return bs
end

# ── Filtered per-cell mutation counts ─────────────────────────────────────────
# Excludes mutations from any branch whose live-descendant count > threshold * N.
# Two-pass: (1) post-order descendant count per node, (2) per-leaf accumulation.

function _count_desc!(node::BinaryNode{NonMarkovCell},
                      desc::Dict{BinaryNode{NonMarkovCell}, Int})
    if isnothing(node.left) && isnothing(node.right)
        desc[node] = 1
        return 1
    end
    count = 0
    isnothing(node.left)  || (count += _count_desc!(node.left,  desc))
    isnothing(node.right) || (count += _count_desc!(node.right, desc))
    desc[node] = count
    return count
end

function compute_filtered_mut_per_cell(root::BinaryNode{NonMarkovCell},
                                       threshold::Float64)
    N         = length(collect(Leaves(root)))
    max_count = floor(Int, threshold * N)

    desc = Dict{BinaryNode{NonMarkovCell}, Int}()
    _count_desc!(root, desc)

    result = Int[]
    for leaf in Leaves(root)
        muts = leaf.data.mutations
        node = leaf
        while !isnothing(node.parent)
            node = node.parent
            desc[node] <= max_count && (muts += node.data.mutations)
        end
        push!(result, muts)
    end
    return result
end
