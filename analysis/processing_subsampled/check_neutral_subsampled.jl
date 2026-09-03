using Pkg
Pkg.activate(dirname(dirname(@__DIR__)))

using MutationLoadDynamics
using Serialization
using Statistics

const HELPERS = joinpath(dirname(@__DIR__), "helpers")
include(joinpath(HELPERS, "types.jl"))

const ROOT = dirname(dirname(@__DIR__))
const FULL = joinpath(ROOT, "data", "processed",            "neutral")
const SUB  = joinpath(ROOT, "data", "processed_subsampled", "neutral")

# Cross-check the subsampled arrays against the full-tree arrays they are meant to
# be co-indexed with. Most of what follows is a property that must hold by
# construction if the two stages are still co-indexed; see the note at the params
# assertion below for the one check that is weaker than it looks.

for (model, stem, N, n) in (
        ("gamma",         "neutral_gamma_N10000_d0.5_k5.0",   10_000, 1_000),
        ("gamma",         "neutral_gamma_N1000_d0.0_k5.0",     1_000,   100),
        ("markov",        "neutral_markov_N10000_d0.9",       10_000,   100),
        ("deterministic", "neutral_deterministic_N16384",     16_384, 1_638),
    )
    fp = q -> joinpath(FULL, model, q, stem * ".jls")
    sp = q -> joinpath(SUB,  model, q, stem * "_n$(n).jls")

    f_sfs, f_mpc, f_ld = deserialize(fp("sfs")), deserialize(fp("mut_per_cell")), deserialize(fp("leaf_depths"))
    s_sfs, s_mpc, s_ld = deserialize(sp("sfs")), deserialize(sp("mut_per_cell")), deserialize(sp("leaf_depths"))
    f_par, s_par = deserialize(fp("params")), deserialize(sp("params"))

    @assert length(s_sfs) == length(f_sfs)                 "sim count differs: $stem"
    # Catches a wrong-shard pairing and a length mismatch. It does NOT catch a
    # reordering for the neutral scenario: every neutral `SimParams` array is 200
    # structurally identical elements (unlike `Sel1Params`/`Sel2Params`, which carry
    # per-simulation `rep`/`seed`/`N_critic`), so a permutation of `s_par` would still
    # equal `f_par` here.
    @assert s_par == f_par                                 "params not co-indexed: $stem"
    @assert all(length.(s_sfs) .== n)                      "sfs length != n: $stem"
    @assert all(length.(s_mpc) .== n)                      "mut_per_cell length != n: $stem"
    @assert all(length.(s_ld)  .== n)                      "leaf_depths length != n: $stem"
    @assert all(length.(f_sfs) .== N)                      "full sfs length != N: $stem"

    # Every mutation counted once per carrier, on both sides of the SFS.
    for i in eachindex(s_sfs)
        @assert sum(k * s_sfs[i][k] for k in 1:n) == sum(s_mpc[i]) "sfs/burden mismatch: $stem sim $i"
    end

    # Sampling can only lose mutations and can never invent a depth or a burden.
    @assert all(sum(s_sfs[i]) <= sum(f_sfs[i]) for i in eachindex(s_sfs)) "sub sfs exceeds full: $stem"
    # Vacuous for the deterministic model: Dirac division timing with no death gives
    # a perfect binary tree, so every leaf sits at the same depth and both sets are
    # the same singleton — this check is only substantive for gamma/markov. The
    # `mut_per_cell` subset check below and the sfs/burden bookkeeping assertion
    # above remain substantive for deterministic too.
    @assert all(issubset(Set(s_ld[i]),  Set(f_ld[i]))  for i in eachindex(s_ld))  "unseen depth: $stem"
    @assert all(issubset(Set(s_mpc[i]), Set(f_mpc[i])) for i in eachindex(s_mpc)) "unseen burden: $stem"

    println(rpad("$model/$stem n=$n", 52),
            " OK   mean burden full=", round(mean(mean.(f_mpc)); digits = 2),
            " sub=", round(mean(mean.(s_mpc)); digits = 2),
            "   segregating sites full=", round(mean(sum.(f_sfs)); digits = 1),
            " sub=", round(mean(sum.(s_sfs)); digits = 1))
end

println("\nAll cross-checks passed.")
