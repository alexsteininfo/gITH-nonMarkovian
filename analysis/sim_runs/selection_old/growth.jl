using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Distributions
using Serialization
using Random
using Statistics

# ── Parameters ────────────────────────────────────────────────────────────────
# Additive selection model: each driver mutation adds a fixed fitness increment δ = s.
# New fitness: f ← f + s. Birth time ~ Gamma(k, 1/(k·f)), so faster division
# is proportional to accumulated fitness. Equivalent to the classic birth-rate
# model b·(1 + accumulated_s) with Gamma timing noise replacing exponential.

const N_TARGET  = 10_000
const N_SIMS    = 50
const b         = 1.0    # baseline birth rate (mean division time = 1/b at f=1)
const d         = 0.5    # death rate
const k         = 5.0    # Gamma shape; CV = 1/√k ≈ 0.45
const ν         = 0.2    # Poisson mean driver mutations per daughter cell
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
        driver_dist           = Dirac(s),           # fixed fitness increment
        fitness_update        = (f, δ) -> f + δ,   # additive accumulation
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
    println("Running $N_SIMS simulations for s = $s (additive, fixed δ = s) ...")

    results = Vector{GrowthSimResult}(undef, N_SIMS)
    Threads.@threads for i in 1:N_SIMS
        results[i] = run_simulation(s, MersenneTwister(i))
    end

    outfile = joinpath(@__DIR__, "..", "..", "data", "raw", "growth",
                       "growth_s$(s)_k$(k).jls")
    mkpath(dirname(outfile))
    serialize(outfile, results)
    println("  → saved $(length(results)) simulations to $outfile")
end

println("All done.")
