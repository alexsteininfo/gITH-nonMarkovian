using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using AbstractTrees
using Serialization
using Random

include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "subsampling.jl"))

const RAW = joinpath(DATA, "raw",            "selection_1")
const OUT = joinpath(DATA, "raw_subsampled", "selection_1")

# Stage 1b for selection scenario 1. The grids mirror
# `analysis/processing/process_sel1.jl` exactly, including `collect(0.1:0.1:2.0)`
# rather than a literal list, so the interpolated filenames match the 280 shards on
# disk character for character.

const S_VALUES = collect(0.1:0.1:2.0)

println("\n── Gamma model ────────────────────────────────────────────────────────")
for N in (1_000, 10_000), d in (0.0, 0.5, 0.9), s in S_VALUES
    raw_path = joinpath(RAW, "gamma", "sel1_gamma_N$(N)_d$(d)_k5.0_s$(s).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    subsample_file(raw_path, joinpath(OUT, "gamma"), sample_sizes(N))
end

println("\n── Markov model ───────────────────────────────────────────────────────")
for N in (1_000, 10_000), d in (0.0, 0.5, 0.9), s in S_VALUES
    raw_path = joinpath(RAW, "markov", "sel1_markov_N$(N)_d$(d)_s$(s).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    subsample_file(raw_path, joinpath(OUT, "markov"), sample_sizes(N))
end

println("\n── Deterministic model ────────────────────────────────────────────────")
for N in (1_024, 16_384), s in S_VALUES
    raw_path = joinpath(RAW, "deterministic", "sel1_deterministic_N$(N)_s$(s).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    subsample_file(raw_path, joinpath(OUT, "deterministic"), sample_sizes(N))
end

println("\nAll done.")
