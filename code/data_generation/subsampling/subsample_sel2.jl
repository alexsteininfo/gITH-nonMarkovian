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

include(joinpath(HELPERS, "types_selection2.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "subsampling.jl"))

const RAW = joinpath(DATA, "raw",            "selection_2")
const OUT = joinpath(DATA, "raw_subsampled", "selection_2")

# Stage 1b for selection scenario 2. The grids mirror
# `analysis/processing/process_sel2.jl` exactly — note that `0.10` and `0.20`
# interpolate as `s0.1` and `s0.2` — so the filenames match the shards on disk
# character for character.

const S_VALUES = [0.05, 0.10, 0.15, 0.20]
const M_CAP    = 10.0

println("\n── Gamma model ────────────────────────────────────────────────────────")
for N in (1_000, 10_000), d in (0.0, 0.5, 0.9), s in S_VALUES
    raw_path = joinpath(RAW, "gamma", "sel2_gamma_N$(N)_d$(d)_k5.0_s$(s)_M$(M_CAP).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    subsample_file(raw_path, joinpath(OUT, "gamma"), sample_sizes(N))
end

println("\n── Markov model ───────────────────────────────────────────────────────")
for N in (1_000, 10_000), d in (0.0, 0.5, 0.9), s in S_VALUES
    raw_path = joinpath(RAW, "markov", "sel2_markov_N$(N)_d$(d)_s$(s)_M$(M_CAP).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    subsample_file(raw_path, joinpath(OUT, "markov"), sample_sizes(N))
end

println("\n── Deterministic model ────────────────────────────────────────────────")
for N in (1_024, 16_384), s in S_VALUES
    raw_path = joinpath(RAW, "deterministic", "sel2_deterministic_N$(N)_s$(s)_M$(M_CAP).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    subsample_file(raw_path, joinpath(OUT, "deterministic"), sample_sizes(N))
end

println("\nAll done.")
