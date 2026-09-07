using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using CopyNumberEvolution

include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "subsampling.jl"))       # for sample_sizes
include(joinpath(HELPERS, "cn_evolution.jl"))

const RAW = joinpath(DATA, "raw_subsampled", "neutral")
const OUT = joinpath(DATA, "CN_subsampled",  "neutral")

# Stage 4 for the neutral scenario: all three timing models — neutral is small
# enough by shard count that none needed to be dropped for this first pass (see
# docs/superpowers/specs/2026-09-07-cn-evolution-design.md). Only the smallest
# sample size of each shard family is used (`last(sample_sizes(N))`).

println("\n── Gamma model ────────────────────────────────────────────────────────")
for N in (1_000, 10_000), d in (0.0, 0.5, 0.9)
    n = last(sample_sizes(N))
    raw_path = joinpath(RAW, "gamma", "neutral_gamma_N$(N)_d$(d)_k5.0_n$(n).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    cn_evolve_file(raw_path, joinpath(OUT, "gamma"))
end

println("\n── Markov model ───────────────────────────────────────────────────────")
for N in (1_000, 10_000), d in (0.0, 0.5, 0.9)
    n = last(sample_sizes(N))
    raw_path = joinpath(RAW, "markov", "neutral_markov_N$(N)_d$(d)_n$(n).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    cn_evolve_file(raw_path, joinpath(OUT, "markov"))
end

println("\n── Deterministic model ────────────────────────────────────────────────")
for N in (1_024, 16_384)
    n = last(sample_sizes(N))
    raw_path = joinpath(RAW, "deterministic", "neutral_deterministic_N$(N)_n$(n).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    cn_evolve_file(raw_path, joinpath(OUT, "deterministic"))
end

println("\nAll done.")
