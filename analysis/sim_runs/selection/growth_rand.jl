using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Distributions
using Serialization
using Random
using Statistics

# ── Parameters ────────────────────────────────────────────────────────────────
# Max-random (clonal competition) selection model.
# Each driver mutation draws a new candidate fitness from Exp(s) + 1.
# The cell adopts the better of its current fitness and the new draw:
#   f ← max(f, 1 + δ),   δ ~ Exp(s)
# Fitness can only increase; once a high-fitness mutation arises, it is locked
# in. Models a scenario where later mutations cannot "undo" earlier gains, and
# only the largest beneficial effect matters.
# For s = 0: Dirac(0) is used so no fitness change occurs (neutral).

const N_TARGET  = 1_000
const N_SIMS    = 50
const b         = 1.0
const d         = 0.5
const k         = 5.0
const ν         = 0.2
const S_VALUES  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

include(joinpath(dirname(@__DIR__), "helpers", "types.jl"))

# ── Per-simulation runner ─────────────────────────────────────────────────────

function run_simulation(s::Float64, rng::AbstractRNG)
    spec = MeasurementSpec(
        trajectory_dt     = 0.1,
        snapshot_triggers = [AtEnd()],
        snapshot_stats    = [],
    )

    # Exponential(0) is undefined; use Dirac(0) for neutral case (no fitness change)
    dist = s > 0.0 ? Exponential(s) : Dirac(0.0)

    pop = initialize_population(fitness_init = 1.0)
    block = NonMarkovBlock(
        birth_dist            = f -> Gamma(k, 1.0 / (k * f)),
        death_dist            = _ -> Gamma(k, 1.0 / (k * d)),
        stopfunction          = pop -> popsize(pop) >= N_TARGET,
        driver_dist           = dist,
        fitness_update        = (f, δ) -> max(f, 1.0 + δ),  # take the best
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
        SimParams(b, d, k, ν, s, N_TARGET, NaN),
    )
end

# ── Main simulation loop ──────────────────────────────────────────────────────

for s in S_VALUES
    println("Running $N_SIMS simulations for s = $s (max-random, δ ~ Exp(s)) ...")

    results = Vector{GrowthSimResult}(undef, N_SIMS)
    Threads.@threads for i in 1:N_SIMS
        results[i] = run_simulation(s, MersenneTwister(i))
    end

    outfile = joinpath(@__DIR__, "..", "..", "data", "raw", "growth_rand",
                       "growth_rand_s$(s)_k$(k).jls")
    mkpath(dirname(outfile))
    serialize(outfile, results)
    println("  → saved $(length(results)) simulations to $outfile")
end

println("All done.")
