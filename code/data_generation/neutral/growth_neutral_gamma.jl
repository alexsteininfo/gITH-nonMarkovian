using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using Distributions
using Serialization
using Random

# ── Parameters ────────────────────────────────────────────────────────────────

const N_SIMS   = 200
const b        = 1.0   # birth rate; mean division time = 1/b
const k        = 5.0   # Gamma shape; CV = 1/√k ≈ 0.45 (cancer-cell regime)
const N_VALUES = [1_000, 10_000]
const D_VALUES = [0.0, 0.5, 0.9]

include(joinpath(HELPERS, "types.jl"))

# ── Per-simulation runner ─────────────────────────────────────────────────────

function run_simulation(N_target::Int, d::Float64, rng::AbstractRNG)
    spec = MeasurementSpec(
        trajectory_dt     = 0.1,
        snapshot_triggers = [AtEnd()],
        snapshot_stats    = [],
    )

    death_scale = d > 0.0 ? 1.0 / (k * d) : 1e8

    pop   = initialize_population(fitness_init = 1.0)
    block = NonMarkovBlock(
        birth_dist            = _ -> Gamma(k, 1.0 / (k * b)),
        death_dist            = _ -> Gamma(k, death_scale),
        stopfunction          = pop -> popsize(pop) >= N_target,
        driver_dist           = Dirac(0.0),
        fitness_update        = (f, δ) -> f,
        ν                     = 2.0,
        restart_on_extinction = true,
    )

    acc = MeasurementAccumulator(spec)
    simulate!(pop, block, rng; accumulator = acc)
    m         = finalize_measurements(acc)
    tree_root = getsingleroot(allcells(pop))

    return GrowthSimResult(
        m.trajectory,
        tree_root,
        SimParams(b, d, k, 2.0, 0.0, N_target, NaN),
    )
end

# ── Main loop ─────────────────────────────────────────────────────────────────

for N_target in N_VALUES, d in D_VALUES
    println("Running $N_SIMS simulations  N=$N_target  d=$d  k=$k ...")

    results = Vector{GrowthSimResult}(undef, N_SIMS)
    Threads.@threads for i in 1:N_SIMS
        results[i] = run_simulation(N_target, d, MersenneTwister(i))
    end

    outfile = joinpath(@__DIR__, "..", "..", "..", "data", "raw", "neutral", "gamma",
                       "neutral_gamma_N$(N_target)_d$(d)_k$(k).jls")
    mkpath(dirname(outfile))
    serialize(outfile, results)
    println("  → $(outfile)")
end

println("All done.")
