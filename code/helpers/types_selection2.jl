# Data structures for selection scenario 2: every mutation carries a fitness effect
# drawn from a Gamma distribution, with a hard cap M on attainable fitness.
#
# Deliberately separate from `types.jl` and `types_selection1.jl`. Julia's
# `Serialization` resolves struct layout at deserialization time, so widening an
# existing struct would make every already-written `.jls` unreadable — 1.8 GB of
# neutral data plus whatever scenario 1 has produced. Each scenario owns its structs.
#
# `include`d (not imported) by every scenario-2 stage so that independent scripts see
# identical definitions. TrajectoryPoint, NonMarkovCell and BinaryNode are exported
# from MutationLoadDynamics.

"""
Design point of one scenario-2 simulation: the background process, the selection
parameters, and the `(model, N_target, d, s, rep)` grid cell it belongs to.

Note the two Gamma shapes are unrelated and both recorded: `gamma_shape` is the shape
of the *division-timing* distribution (`CV = 1/√gamma_shape`), while `effect_shape` is
the shape of the *per-mutation effect* distribution.
"""
struct Sel2Params
    b::Float64            # baseline birth rate at f = 1
    d::Float64            # death rate (0.0 = no death); fitness-independent
    gamma_shape::Float64  # timing shape k: 5.0 gamma, 1.0 markov, Inf deterministic
    nu::Float64           # mutations per daughter per division — all of them drivers
    N_target::Int         # stop condition
    model::Symbol         # :deterministic | :gamma | :markov
    s::Float64            # mean per-mutation selection coefficient
    M::Float64            # fitness cap
    effect_shape::Float64 # a of the effect distribution Gamma(a, 1/a): mean 1, CV 1/√a
    rep::Int              # replicate index within the (model, N_target, d, s) cell
    seed::UInt64          # rng seed, so any single simulation replays in isolation
end

"""
Outcome of one simulation. `n_restarts` counts how many times the population went
extinct before the attempt that reached `N_target`; `restart_on_extinction` discards
the failed attempt's trajectory, so what is stored covers the surviving run only.

`1 / (1 + mean(n_restarts))` over a cell estimates the survival probability of a
single founder — itself a result, since Gamma and exponential timing give different
extinction probabilities at the same `b` and `d`.

Every leaf's fitness already lives in `NonMarkovCell.fitness`, so `tree_root` carries
the full final fitness distribution and no snapshot statistics are needed.
"""
struct Sel2SimResult
    trajectory::Vector{TrajectoryPoint}
    tree_root::Union{BinaryNode{NonMarkovCell}, Nothing}
    params::Sel2Params
    n_restarts::Int
end
