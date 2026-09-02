using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Distributions
using Serialization
using Random
using Statistics

# ── Parameters ────────────────────────────────────────────────────────────────
# Multiplicative selection model with random effect (X ~ Exp(1)).
# Each driver mutation draws X ~ Exp(1) and updates fitness:
#   f ← min( f · (1 + s · X · (1 − f/M)),  M )
# The expected per-mutation gain is E[sX(1−f/M)] = s(1−f/M), the same as
# the fixed-effect model, but with additional variance from X. This models
# driver mutations whose functional magnitude is unknown and exponentially
# distributed with mean 1. For s = 0 the rule reduces to f ← f (neutral).

const N_TARGET  = 1_000
const N_SIMS    = 50
const b         = 1.0
const d         = 0.5
const k         = 5.0
const ν         = 0.2
const M         = 10.0
const S_VALUES  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

include(joinpath(dirname(@__DIR__), "helpers", "types.jl"))

# ── Per-simulation runner ─────────────────────────────────────────────────────

function run_simulation(s::Float64, rng::AbstractRNG)
    spec = MeasurementSpec(
        trajectory_dt     = 0.1,
        snapshot_triggers = [AtEnd()],
        snapshot_stats    = [],
    )

    pop = initialize_population(fitness_init = 1.0)
    block = NonMarkovBlock(
        birth_dist            = f -> Gamma(k, 1.0 / (k * f)),
        death_dist            = _ -> Gamma(k, 1.0 / (k * d)),
        stopfunction          = pop -> popsize(pop) >= N_TARGET,
        driver_dist           = Exponential(1.0),                        # X ~ Exp(1)
        fitness_update        = (f, δ) -> min(f * (1.0 + s * δ * (1.0 - f / M)), M),
        ν                     = ν,
        restart_on_extinction = true,
    )

    acc = MeasurementAccumulator(spec)
    simulate!(pop, block, rng; accumulator = acc)
    m         = finalize_measurements(acc)
    tree_root = getsingleroot(allcells(pop))

    return GrowthSimResult(
        m.trajectory,
        tree_root,
        SimParams(b, d, k, ν, s, N_TARGET, M),
    )
end

# ── Main simulation loop ──────────────────────────────────────────────────────

for s in S_VALUES
    println("Running $N_SIMS simulations for s = $s (multiplicative, X ~ Exp(1)) ...")

    results = Vector{GrowthSimResult}(undef, N_SIMS)
    Threads.@threads for i in 1:N_SIMS
        results[i] = run_simulation(s, MersenneTwister(i))
    end

    outfile = joinpath(@__DIR__, "..", "..", "data", "raw", "growth_multi_rand",
                       "growth_multi_rand_s$(s)_k$(k).jls")
    mkpath(dirname(outfile))
    serialize(outfile, results)
    println("  → saved $(length(results)) simulations to $outfile")
end

println("All done.")
