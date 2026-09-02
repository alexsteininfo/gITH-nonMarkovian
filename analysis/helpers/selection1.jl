# Shared machinery for selection scenario 1. Lives here rather than in the three
# `analysis/sim_runs/selection_1/growth_sel1_*.jl` scripts so the accept/retry loop
# exists once instead of three times; those scripts supply only the model-specific
# birth/death distributions and the output path.
#
# The calling script must already have done `using MutationLoadDynamics,
# Distributions, Random, Statistics` and included `types_selection1.jl`.

"""
    ncritic_grid(N_target) -> Vector{Int}

Ten log-spaced injection sizes from 1 to `N_target ÷ 2`.

Log spacing makes the injection *times* approximately equally spaced: during
exponential growth `N(t) ≈ e^{rt}`, so `t ≈ ln(N)/r`. Ten values rather than
twenty because twenty log-spaced integers over `[1, 500]` collide at the low end
(1, 1, 2, 3, 4, 5, …) and dedup down to about seventeen.
"""
ncritic_grid(N_target::Int) =
    unique(round.(Int, exp.(range(0.0, log(N_target ÷ 2), length = 10))))

"""
    run_sel1_once(birth_dist, death_dist, N_target, s, N_critic, rng;
                  nu = 2.0, trajectory_dt = 0.1)

Run one simulation to `N_target`, injecting a single driver of strength `s` into the
left daughter of the division that first takes the population from `N_critic` to
`N_critic + 1` alive cells.

The injection uses the `on_division` hook (MutationLoadDynamics ≥ 0.2.0), which fires
after `celldivision!` and **before** either daughter is scheduled — so the boosted
daughter's own first division is drawn at `f = 1 + s`, not just its descendants'.
Because `celldivision!` removes the parent and inserts both daughters before the hook
runs, the trigger inside the callback is `popsize(pop) == N_critic + 1`.

Makes no judgement about whether the run is usable; that is `run_sel1_accepted`'s job.
Returns a NamedTuple: `pop, trajectory, injected, t_inject, N_at_inject,
driver_cell_id, driver_node`.
"""
function run_sel1_once(
    birth_dist,
    death_dist,
    N_target::Int,
    s::Float64,
    N_critic::Int,
    rng::AbstractRNG;
    nu::Float64 = 2.0,
    trajectory_dt::Float64 = 0.1,
)
    injected    = Ref(false)
    t_inject    = Ref(NaN)
    N_at_inject = Ref(0)
    driver_id   = Ref(zero(Int64))
    driver_node = Ref{Union{Nothing, BinaryNode{NonMarkovCell}}}(nothing)

    spec = MeasurementSpec(
        trajectory_dt     = trajectory_dt,
        snapshot_triggers = [AtEnd()],
        snapshot_stats    = [],
    )

    pop   = initialize_population(fitness_init = 1.0)
    block = NonMarkovBlock(
        birth_dist            = birth_dist,
        death_dist            = death_dist,
        stopfunction          = pop -> popsize(pop) >= N_target,
        driver_dist           = Dirac(0.0),   # neutral channel: no fitness effect
        fitness_update        = (f, δ) -> f,  # neutral channel: fitness untouched
        ν                     = nu,
        restart_on_extinction = false,        # the retry loop owns extinction
        on_division = function (pop, parent, d1, d2)
            injected[] && return nothing
            popsize(pop) == N_critic + 1 || return nothing
            injected[] = true
            c = d1.data
            # NonMarkovCell is immutable, BinaryNode is not: replace .data wholesale.
            # `mutations` is carried over unchanged so the driver stays out of the
            # neutral SFS channel; the clone is identified by fitness > 1.
            d1.data = NonMarkovCell(c.id, c.birthtime, c.mutations, 1.0 + s)
            t_inject[]    = c.birthtime
            N_at_inject[] = popsize(pop)
            driver_id[]   = c.id
            driver_node[] = d1
            return nothing
        end,
    )

    acc = MeasurementAccumulator(spec)
    simulate!(pop, block, rng; accumulator = acc)
    m = finalize_measurements(acc)

    return (
        pop            = pop,
        trajectory     = m.trajectory,
        injected       = injected[],
        t_inject       = t_inject[],
        N_at_inject    = N_at_inject[],
        driver_cell_id = driver_id[],
        driver_node    = driver_node[],
    )
end
