using Pkg
const ROOT = dirname(dirname(dirname(dirname(@__DIR__))))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using Distributions
using Serialization
using Statistics
using Random

# ── Parameters ────────────────────────────────────────────────────────────────
# Exponential waiting times (Gamma with k=1): the Markovian / Gillespie limit.
# Background process identical to analysis/sim_runs/neutral/growth_neutral_markov.jl.

const b        = 1.0
const k        = 1.0   # recorded as gamma_shape; exponential is the k=1 Gamma
const nu       = 2.0
const N_VALUES = [1_000, 10_000]
const D_VALUES = [0.0, 0.5, 0.9]
const S_VALUES = collect(0.1:0.1:2.0)
const N_REPS   = 5

include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "selection1.jl"))

# ── Main loop ─────────────────────────────────────────────────────────────────

for N_target in N_VALUES, d in D_VALUES
    birth_dist = f -> Exponential(1.0 / (b * f))
    death_dist = _ -> Exponential(d > 0.0 ? 1.0 / d : 1e8)
    N_CRITIC   = ncritic_grid(N_target)

    for s in S_VALUES
        # The sweep is resumable: a shard is skipped once its output file exists, so a
        # crash mid-sweep restarts only from the first missing shard. Delete a shard's
        # .jls file to force it to be regenerated.
        outfile = joinpath(@__DIR__, "..", "..", "..", "data", "raw", "selection_1",
                           "markov", "sel1_markov_N$(N_target)_d$(d)_s$(s).jls")
        if isfile(outfile)
            println("  → $(outfile) exists, skipping")
            continue
        end

        println("markov  N=$N_target  d=$d  s=$s ...")

        cells   = [(nc, rep) for nc in N_CRITIC for rep in 1:N_REPS]
        results = Vector{Sel1SimResult}(undef, length(cells))
        Threads.@threads for i in eachindex(cells)
            nc, rep = cells[i]
            results[i] = run_sel1_accepted(birth_dist, death_dist,
                Sel1Params(b, d, k, nu, N_target, :markov, s, nc, rep))
        end

        mkpath(dirname(outfile))
        serialize(outfile, results)
        println("  → $(outfile)  (mean attempts $(round(mean(r.injection.n_attempts for r in results), digits = 2)))")
    end
end

println("All done.")
