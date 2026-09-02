# Data structures for selection scenario 1: a single fitness-enhancing mutation
# injected at a prescribed population size.
#
# Deliberately separate from `types.jl`. Julia's `Serialization` resolves struct
# layout at deserialization time, so adding a field to `SimParams` or
# `GrowthSimResult` would make every already-written neutral `.jls` unreadable.
# Scenario 1 therefore gets its own structs rather than widening the neutral ones.
#
# `include`d (not imported) by every scenario-1 stage so that independent scripts
# see identical definitions. TrajectoryPoint, NonMarkovCell and BinaryNode are
# exported from MutationLoadDynamics.

"""
Design point of one scenario-1 simulation: the background process plus the
`(s, N_critic, rep)` grid cell it belongs to.
"""
struct Sel1Params
    b::Float64            # baseline birth rate at f = 1
    d::Float64            # death rate (0.0 = no death)
    gamma_shape::Float64  # k: 5.0 gamma, 1.0 markov, Inf deterministic
    nu::Float64           # neutral mutations per daughter per division (2.0)
    N_target::Int         # stop condition
    model::Symbol         # :deterministic | :gamma | :markov
    s::Float64            # selection coefficient of the single driver
    N_critic::Int         # popsize before the injecting division
    rep::Int              # replicate index within the (s, N_critic) cell
end

"""
Measured outcome of the injection. `N_at_inject` is recorded rather than assumed
equal to `N_critic + 1`, so a future change of mechanism is detectable from the
data alone.
"""
struct Sel1Injection
    t_inject::Float64      # simulation time of the injecting division
    N_at_inject::Int       # popsize immediately after it
    driver_cell_id::Int64  # id of the boosted daughter
    driver_clone_size::Int # alive cells with fitness > 1 at N_target
    n_attempts::Int        # attempts needed before acceptance
    seed::UInt64           # rng seed of the accepted attempt
end

struct Sel1SimResult
    trajectory::Vector{TrajectoryPoint}
    tree_root::Union{BinaryNode{NonMarkovCell}, Nothing}
    params::Sel1Params
    injection::Sel1Injection
end
