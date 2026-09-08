using Pkg
const ROOT = dirname(dirname(dirname(dirname(@__DIR__))))
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
const b        = 1.0   # birth rate; every cell divides at exactly time 1/b
const N_VALUES = [1_024, 16_384]

include(joinpath(HELPERS, "types.jl"))

# ── Per-simulation runner ─────────────────────────────────────────────────────

function run_simulation(N_target::Int, rng::AbstractRNG)
    spec = MeasurementSpec(
        trajectory_dt     = 0.1,
        snapshot_triggers = [AtEnd()],
        snapshot_stats    = [],
    )

    # Dirac waiting times: every cell divides at exactly 1/b after birth.
    # No death (d=0): the resulting tree is a perfect binary tree.
    # Population size overshoots to the next power of 2 above N_target.
    pop   = initialize_population(fitness_init = 1.0)
    block = NonMarkovBlock(
        birth_dist            = _ -> Dirac(1.0 / b),
        death_dist            = _ -> Dirac(1e8),
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
        SimParams(b, 0.0, Inf, 2.0, 0.0, N_target, NaN),
    )
end

# ── Main loop ─────────────────────────────────────────────────────────────────

for N_target in N_VALUES
    println("Running $N_SIMS simulations  N=$N_target  (deterministic / Dirac) ...")

    results = Vector{GrowthSimResult}(undef, N_SIMS)
    Threads.@threads for i in 1:N_SIMS
        results[i] = run_simulation(N_target, MersenneTwister(i))
    end

    outfile = joinpath(DATA, "raw", "neutral", "deterministic",
                       "neutral_deterministic_N$(N_target).jls")
    mkpath(dirname(outfile))
    serialize(outfile, results)
    println("  → $(outfile)")
end

println("All done.")
