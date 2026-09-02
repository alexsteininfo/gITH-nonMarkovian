using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using AbstractTrees
using Distributions
using Serialization
using Statistics
using Random
using Test

const ROOT    = dirname(dirname(dirname(@__DIR__)))
const HELPERS = joinpath(dirname(dirname(@__DIR__)), "helpers")

include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_selection1.jl"))

@testset "selection scenario 1" begin

@testset "types round-trip through Serialization" begin
    p   = Sel1Params(1.0, 0.5, 5.0, 2.0, 1000, :gamma, 0.3, 44, 2)
    inj = Sel1Injection(3.5, 45, 87, 120, 4, 0x1234)
    r   = Sel1SimResult(TrajectoryPoint[], nothing, p, inj)

    path = tempname() * ".jls"
    serialize(path, [r])
    back = deserialize(path)::Vector{Sel1SimResult}
    @test length(back) == 1
    @test back[1].params    == p
    @test back[1].injection == inj
    @test back[1].params.model    === :gamma
    @test back[1].params.N_critic == 44
    rm(path)
end

@testset "existing neutral .jls remain readable" begin
    f = joinpath(ROOT, "data", "raw", "neutral", "gamma",
                 "neutral_gamma_N1000_d0.5_k5.0.jls")
    if isfile(f)
        sims = deserialize(f)
        @test sims isa Vector{GrowthSimResult}
        @test sims[1].params isa SimParams
        @test sims[1].params.nu == 2.0
    else
        @info "skipping: $f not present"
    end
end

end
