# Shared data structures for growth simulations.
# Included by both sim scripts and plot scripts so that Julia's Serialization
# module can resolve the types in both contexts.
# TrajectoryPoint, NonMarkovCell, BinaryNode are exported from MutationLoadDynamics.

struct SimParams
    b::Float64            # background birth rate (mean division time = 1/b)
    d::Float64            # death rate (mean death time = 1/d; 0.0 = no death)
    gamma_shape::Float64  # k for Gamma(k, θ) distributions (CV = 1/√k)
    nu::Float64           # driver mutation rate per daughter cell (0.0 = neutral)
    s::Float64            # selection coefficient for this run (0.0 = neutral)
    N_target::Int         # target final population size
    M::Float64            # max birth rate cap for multiplicative models (NaN = not applicable)
end

struct GrowthSimResult
    trajectory::Vector{TrajectoryPoint}
    tree_root::Union{BinaryNode{NonMarkovCell}, Nothing}
    params::SimParams
end
