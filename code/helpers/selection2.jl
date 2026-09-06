# Shared machinery for selection scenario 2. Lives here rather than in the three
# `analysis/sim_runs/selection_2/growth_sel2_*.jl` scripts so the effect distribution,
# the fitness-update rule and the seeding exist once instead of three times; those
# scripts supply only the model-specific birth/death distributions and the output path.
#
# The calling script must already have done `using MutationLoadDynamics, Distributions,
# Random` and included `types_selection2.jl`.
#
# Scenario 2 vs. scenario 1: there is no injection hook and no accept/retry loop.
# Every mutation is a driver, so no single mutation's fate can make a run unusable —
# subclones are lost to drift constantly and that is the process, not a failure.
# `restart_on_extinction = true` therefore suffices, exactly as in the neutral runs;
# it fires only on *total* population extinction (`N == 0`).

"""
    sel2_effect_dist(a) -> Gamma

Distribution of the per-mutation effect multiplier `X`, normalised to `E[X] = 1` so
that `s` alone sets the mean per-mutation selection coefficient and `a` sets only its
spread (`CV = 1/√a`). `a = 1` recovers the `Exponential(1)` of the superseded
`selection_old/growth_multi_rand.jl`.
"""
sel2_effect_dist(a::Float64) = Gamma(a, 1.0 / a)

"""
    sel2_fitness_update(s, M) -> (f, X) -> Float64

The multiplicative, capped update applied once per mutation:

    f ← min( f · (1 + s · X · (1 − f/M)),  M )

Same functional form as `selection_old/growth_multi_rand.jl`. The logistic factor
`(1 − f/M)` shrinks the gain as fitness approaches the cap `M`, so per-mutation
returns diminish rather than compounding without bound; the outer `min` keeps `f ≤ M`
even for a large draw of `X`. Fitness is non-decreasing — there are no deleterious
mutations in this scenario. At `s = 0` the rule reduces to `f ← f` exactly.

Because `E[X] = 1`, the expected gain is `s·f·(1 − f/M)` regardless of `effect_shape`,
which is what makes an `effect_shape` sweep a pure variance comparison at fixed mean.
"""
sel2_fitness_update(s::Float64, M::Float64) =
    (f, X) -> min(f * (1.0 + s * X * (1.0 - f / M)), M)

"""
    sel2_seed(model, N_target, d, s, rep) -> UInt64

Deterministic per-simulation seed. Derived from the whole grid cell rather than from
`rep` alone so that different cells do not share rng streams, and recorded on the
result so any single simulation can be reproduced in isolation without replaying the
sweep.
"""
sel2_seed(model::Symbol, N_target::Int, d::Float64, s::Float64, rep::Int) =
    hash((:sel2, model, N_target, d, s, rep))

"""
    run_sel2_once(birth_dist, death_dist, params; trajectory_dt = 0.02) -> Sel2SimResult

Grow one population from a single founder to `params.N_target`, restarting on total
extinction, and package the result.

`birth_dist` and `death_dist` are the model-specific `fitness -> Distribution`
closures. Birth is fitness-scaled in all three timing models; death is not.

`trajectory_dt` defaults to 0.02 rather than the neutral runs' 0.1: with fitness
reaching several times the baseline the population hits `N_target` far sooner in
simulation time, and 0.1 would leave only a couple of dozen trajectory points.

The returned `tree_root` may be `nothing` if `getsingleroot` finds no unique root;
that mirrors the neutral runs, and `process_*` skips such entries.
"""
function run_sel2_once(
    birth_dist,
    death_dist,
    params::Sel2Params;
    trajectory_dt::Float64 = 0.02,
)
    # The package cannot reset closure state on restart, which is exactly what is
    # wanted here: the counter accumulates across attempts.
    n_restarts = Ref(0)

    spec = MeasurementSpec(
        trajectory_dt     = trajectory_dt,
        snapshot_triggers = [AtEnd()],
        snapshot_stats    = [],
    )

    pop   = initialize_population(fitness_init = 1.0)
    block = NonMarkovBlock(
        birth_dist            = birth_dist,
        death_dist            = death_dist,
        stopfunction          = pop -> popsize(pop) >= params.N_target,
        driver_dist           = sel2_effect_dist(params.effect_shape),
        fitness_update        = sel2_fitness_update(params.s, params.M),
        ν                     = params.nu,
        restart_on_extinction = true,
        on_restart            = function (_)
            n_restarts[] += 1
            return nothing
        end,
    )

    acc = MeasurementAccumulator(spec)
    simulate!(pop, block, MersenneTwister(params.seed); accumulator = acc)
    m = finalize_measurements(acc)

    return Sel2SimResult(
        m.trajectory,
        getsingleroot(allcells(pop)),
        params,
        n_restarts[],
    )
end
