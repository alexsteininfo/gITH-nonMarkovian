const ROOT = normpath(joinpath(@__DIR__))
using Pkg; Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using Serialization
using MutationLoadDynamics
using EvoTracer
include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "evotracer_io.jl"))
include(joinpath(HELPERS, "inference_runners.jl"))
include(joinpath(HELPERS, "theory.jl"))

in_path = joinpath(DATA, "raw", "selection_1", "deterministic", "sel1_deterministic_N1024_s0.1.jls")
sims = deserialize(in_path)
sim = sims[1]
if !isnothing(sim.tree_root)
    fresh = selection_suite(prepared_full(sim; primary = :mutations))
    println("Keys in fresh selection_1 output:")
    for k in sort(collect(keys(fresh)))
        println("  :$k")
    end
end
