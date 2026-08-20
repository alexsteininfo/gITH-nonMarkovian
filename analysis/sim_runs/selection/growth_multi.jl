using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Distributions
using Serialization
using Random
using Statistics

# ── Parameters ────────────────────────────────────────────────────────────────
# Multiplicative selection model with fixed effect (X = 1).
# Each driver mutation updates fitness via a logistic-growth rule with cap M:
#   f ← min( f · (1 + s · (1 − f/M)),  M )
# This is the discrete-time analogue of logistic growth toward carrying capacity M.
# The factor (1 − f/M) reduces the per-mutation gain as fitness approaches M,
# preventing unbounded acceleration. For s = 0 the rule reduces to f ← f (neutral).

const N_TARGET  = 1_000
const N_SIMS    = 50
const b         = 1.0
const d         = 0.5
const k         = 5.0
const ν         = 0.2
const M         = 10.0   # maximum birth rate (fitness cap)
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
        driver_dist           = Dirac(1.0),                              # X = 1 (fixed)
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
    println("Running $N_SIMS simulations for s = $s (multiplicative, fixed X = 1) ...")

    results = Vector{GrowthSimResult}(undef, N_SIMS)
    Threads.@threads for i in 1:N_SIMS
        results[i] = run_simulation(s, MersenneTwister(i))
    end

    outfile = joinpath(@__DIR__, "..", "..", "data", "raw", "growth_multi",
                       "growth_multi_s$(s)_k$(k).jls")
    mkpath(dirname(outfile))
    serialize(outfile, results)
    println("  → saved $(length(results)) simulations to $outfile")
end

println("All done.")
