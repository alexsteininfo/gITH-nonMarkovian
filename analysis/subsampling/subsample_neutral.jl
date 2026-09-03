using Pkg
Pkg.activate(dirname(dirname(@__DIR__)))

using MutationLoadDynamics
using AbstractTrees
using Serialization
using Random

const HELPERS = joinpath(dirname(@__DIR__), "helpers")

include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "subsampling.jl"))

# ── Paths ─────────────────────────────────────────────────────────────────────

const RAW = joinpath(@__DIR__, "..", "..", "data", "raw",            "neutral")
const OUT = joinpath(@__DIR__, "..", "..", "data", "raw_subsampled", "neutral")

# ── Subsampling ───────────────────────────────────────────────────────────────
#
# Stage 1b: from each full tree of N cells, draw n cells uniformly without
# replacement and store the induced lineage tree — the sampled leaves plus every
# ancestor of a sampled leaf, unary nodes retained. See `helpers/subsampling.jl`
# for why nothing is collapsed, and the sample sizes in `sample_sizes`.
#
# The grids below mirror `analysis/processing/process_neutral.jl` exactly, so the
# interpolated filenames match the shards on disk character for character.
# `data/raw/neutral/deterministic/` also holds superseded `N1000`/`N10000` files;
# like `process_neutral.jl`, this script ignores them in favour of `N1024`/`N16384`.

# ── Gamma shards (6): N ∈ {1000, 10000}, d ∈ {0.0, 0.5, 0.9}, k = 5.0 ─────────

println("\n── Gamma model ────────────────────────────────────────────────────────")
for N in (1_000, 10_000), d in (0.0, 0.5, 0.9)
    raw_path = joinpath(RAW, "gamma", "neutral_gamma_N$(N)_d$(d)_k5.0.jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    subsample_file(raw_path, joinpath(OUT, "gamma"), sample_sizes(N))
end

# ── Markov shards (6): N ∈ {1000, 10000}, d ∈ {0.0, 0.5, 0.9} ────────────────

println("\n── Markov model ───────────────────────────────────────────────────────")
for N in (1_000, 10_000), d in (0.0, 0.5, 0.9)
    raw_path = joinpath(RAW, "markov", "neutral_markov_N$(N)_d$(d).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    subsample_file(raw_path, joinpath(OUT, "markov"), sample_sizes(N))
end

# ── Deterministic shards (2): N ∈ {1024, 16384}, d = 0.0 ─────────────────────

println("\n── Deterministic model ────────────────────────────────────────────────")
for N in (1_024, 16_384)
    raw_path = joinpath(RAW, "deterministic", "neutral_deterministic_N$(N).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    subsample_file(raw_path, joinpath(OUT, "deterministic"), sample_sizes(N))
end

println("\nAll done.")
