using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using CopyNumberEvolution

include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "subsampling.jl"))       # for sample_sizes
include(joinpath(HELPERS, "cn_evolution.jl"))

const RAW = joinpath(DATA, "raw_subsampled", "selection_1")
const OUT = joinpath(DATA, "CN_subsampled",  "selection_1")

# Stage 4 for selection scenario 1: gamma model only for this first pass (markov
# and deterministic deferred — see the design doc). sel1's gamma grid only ever
# uses N ∈ {1000, 10000}, so every shard here is at n = last(sample_sizes(N)) = 100.

const S_VALUES = collect(0.1:0.1:2.0)

println("\n── Gamma model ────────────────────────────────────────────────────────")
for N in (1_000, 10_000), d in (0.0, 0.5, 0.9), s in S_VALUES
    n = last(sample_sizes(N))
    raw_path = joinpath(RAW, "gamma", "sel1_gamma_N$(N)_d$(d)_k5.0_s$(s)_n$(n).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    cn_evolve_file(raw_path, joinpath(OUT, "gamma"))
end

println("\nAll done.")
