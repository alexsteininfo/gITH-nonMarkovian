# Newick serialization for MutationLoadDynamics.jl lineage trees.
# Include via: include(joinpath(HELPERS, "newick.jl"))
# Requires: using MutationLoadDynamics, AbstractTrees in the calling script.
#
# Deliberately knows nothing about GrowthSimResult / SubsampleResult: it takes a
# BinaryNode and a tmax, so the same writer serves full and subsampled trees.

using Printf

# ── Atomic write ──────────────────────────────────────────────────────────────

"""
    write_atomic(path, content)

Write `content` to `path` via a temporary file plus rename.

The export stage decides "already done" from `isfile`, so a kill mid-write would
otherwise leave a truncated-but-present file that a later run silently reuses. Same
reasoning, and same shape, as `serialize_atomic` in `subsampling.jl`.
"""
function write_atomic(path::AbstractString, content::AbstractString)
    tmp = path * ".tmp"
    open(tmp, "w") do io
        write(io, content)
    end
    mv(tmp, path; force = true)
    return nothing
end

# ── Time span ─────────────────────────────────────────────────────────────────

"""
    tree_tmax(root) -> Float64

The latest leaf birthtime in the tree.

For a full tree this is the exact simulation stop time: the run halts on the division
that pushes the population to `N_target`, so the last-born cells are leaves. It is a
better estimate than `trajectory[end].t`, which is sampled at `dt = 0.1` and therefore
lags (39.000 vs 39.076 in `neutral_markov_N1000_d0.9.jls`). Callers holding a
`GrowthSimResult` should take the max of the two.

For a *subsampled* tree it slightly underestimates the population's stop time, by the
gap between the last birth among `n` sampled cells and among all `N` — of order
`1/(n·ln N)` relative, so ~0.1% at `n = 100`. `verify_tree_export.jl` checks that bound.
"""
tree_tmax(root::BinaryNode{NonMarkovCell}) =
    maximum(l.data.birthtime for l in Leaves(root))

# ── Newick ────────────────────────────────────────────────────────────────────

# 12 significant digits: a root-to-leaf path can be ~100 edges long, and ape's
# is.ultrametric tolerance is ~1.5e-8, so per-edge rounding must stay well below that.
_nwk_len(x::Real) = @sprintf("%.12g", x)

"""
    newick_string(root; tmax) -> String

Newick representation of the lineage tree rooted at `root`.

Branch lengths are in simulation time: the edge *into* a node is that cell's lifetime,
`child.birthtime - node.birthtime`, or `tmax - node.birthtime` for a leaf. Tip labels
are `NonMarkovCell.id`. The root's own lifetime is emitted as a root branch length,
which `ape::read.tree` stores in `tree\$root.edge`; the R side adds it back so tips
land at exactly `t/tmax == 1`.

Because every leaf is alive at `tmax`, the result is **exactly ultrametric** — every
root-to-leaf branch-length sum equals `tmax`. That is the strongest single check on
this writer and both verifiers assert it.

**Unary chains are collapsed.** A node with one surviving child contributes its
lifetime to that child's edge instead of becoming a node of its own. Those nodes are
divisions whose sibling lineage went extinct and was removed by `prune_tree!`; at
`d = 0.9` they are 64-76% of all internal nodes. This is a deliberate local exception
to the "prune, do not collapse" invariant in `CLAUDE.md`: that invariant protects the
stored `.jls` trees, where collapsing would corrupt `leaf_depths` and
`mutations_per_cell`. Here it happens downstream of every statistic, in a file used
only for drawing; it preserves leaf times exactly, shrinks the markov `d = 0.9` tree
from 5197 nodes to 1999, and avoids relying on ape's handling of Newick singletons.
Leaf depth must therefore be measured on the tree, not on the Newick.
"""
function newick_string(root::BinaryNode{NonMarkovCell}; tmax::Real)
    io = IOBuffer()
    _newick!(io, root, Float64(tmax), 0.0)
    print(io, ';')
    return String(take!(io))
end

# `carry` is branch length accumulated from collapsed unary ancestors of `node`.
function _newick!(io::IOBuffer, node::BinaryNode{NonMarkovCell}, tmax::Float64, carry::Float64)
    l, r = node.left, node.right

    if l === nothing && r === nothing
        print(io, node.data.id, ':', _nwk_len(carry + tmax - node.data.birthtime))
        return
    end

    if l === nothing || r === nothing
        child = l === nothing ? r : l
        _newick!(io, child, tmax,
                 carry + child.data.birthtime - node.data.birthtime)
        return
    end

    # Both daughters are created by the same division, so they share a birthtime.
    print(io, '(')
    _newick!(io, l, tmax, 0.0)
    print(io, ',')
    _newick!(io, r, tmax, 0.0)
    print(io, ')', ':', _nwk_len(carry + l.data.birthtime - node.data.birthtime))
    return
end

"""
    write_newick(path, root; tmax)

`newick_string` written atomically to `path`.
"""
write_newick(path::AbstractString, root::BinaryNode{NonMarkovCell}; tmax::Real) =
    write_atomic(path, newick_string(root; tmax = tmax))
