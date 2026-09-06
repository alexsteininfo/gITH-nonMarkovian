using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using AbstractTrees
using Distributions
using Serialization
using Statistics
using Random
using Test

include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "selection1.jl"))

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

@testset "ncritic_grid matches the spec table" begin
    @test ncritic_grid(1_000)  == [1, 2, 4, 8, 16, 32, 63, 126, 251, 500]
    @test ncritic_grid(10_000) == [1, 3, 7, 17, 44, 113, 292, 753, 1941, 5000]
    @test ncritic_grid(1_024)  == [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
    @test ncritic_grid(16_384) == [1, 3, 7, 20, 55, 149, 406, 1106, 3010, 8192]

    for N in (1_000, 10_000, 1_024, 16_384)
        g = ncritic_grid(N)
        @test length(g) == 10          # no collisions after dedup
        @test first(g) == 1
        @test last(g)  == N ÷ 2
        @test issorted(g)
        @test allunique(g)
    end
end

@testset "driver injection mechanics" begin
    birth = f -> Gamma(5.0, 1.0 / (5.0 * f))
    death = _ -> Gamma(5.0, 1.0 / (5.0 * 0.5))

    r = run_sel1_once(birth, death, 1_000, 1.0, 100, MersenneTwister(1))
    @test r.injected
    @test r.N_at_inject == 101            # hook sees popsize == N_critic + 1
    @test popsize(r.pop) >= 1_000
    @test r.driver_cell_id > 0

    # The boosted cells are exactly the alive leaves under the injected node —
    # one clade, no leakage into unrelated lineages.
    clone = count(>(1.0), fitness_per_cell(r.pop))
    @test clone > 0
    @test length(collect(Leaves(r.driver_node))) == clone
end

@testset "boosted cell divides at the boosted rate" begin
    # With Dirac timing a cell of fitness f has lifetime exactly 1/f. This is the
    # test that catches an injection applied after the daughter was scheduled:
    # that ordering would leave the driver's OWN first division at 1/1 = 1.0 and
    # only speed up its descendants.
    s = 1.0
    r = run_sel1_once(f -> Dirac(1.0 / f), _ -> Dirac(1e8), 8, s, 1,
                      MersenneTwister(2); nu = 0.0)
    @test r.injected
    @test celllifetime(r.driver_node) ≈ 1.0 / (1.0 + s)
end

@testset "s = 0 injection does not perturb the rng stream" begin
    # N_critic = 0 can never fire, because popsize inside the hook is always >= 2.
    # So this compares an s=0 injection against an identical un-hooked run.
    birth = f -> Gamma(5.0, 1.0 / (5.0 * f))
    death = _ -> Gamma(5.0, 1.0 / (5.0 * 0.5))
    a = run_sel1_once(birth, death, 500, 0.0, 50, MersenneTwister(99))
    b = run_sel1_once(birth, death, 500, 0.0,  0, MersenneTwister(99))

    @test a.injected
    @test !b.injected
    @test a.pop.t === b.pop.t
    @test popsize(a.pop) == popsize(b.pop)
    @test sort(mutations_per_cell(a.pop)) == sort(mutations_per_cell(b.pop))
    @test all(f === 1.0 for f in fitness_per_cell(a.pop))
end

@testset "acceptance and retry" begin
    birth   = f -> Gamma(5.0, 1.0 / (5.0 * f))
    nodeath = _ -> Gamma(5.0, 1e8)
    death9  = _ -> Gamma(5.0, 1.0 / (5.0 * 0.9))

    # d = 0: no extinction, no driver loss, so the first attempt is always accepted.
    p0 = Sel1Params(1.0, 0.0, 5.0, 2.0, 500, :gamma, 0.5, 50, 1)
    r0 = run_sel1_accepted(birth, nodeath, p0)
    @test r0 isa Sel1SimResult
    @test r0.params == p0
    @test r0.injection.n_attempts == 1
    @test r0.injection.N_at_inject == 51
    @test r0.injection.driver_clone_size > 0
    @test !isnothing(r0.tree_root)
    @test length(r0.trajectory) > 0

    # d = 0.9: extinction and driver loss are both common, so retries must happen.
    # Asserted over ten replicates rather than one, so the check is not a coin flip.
    attempts = [run_sel1_accepted(birth, death9,
                    Sel1Params(1.0, 0.9, 5.0, 2.0, 1_000, :gamma, 0.1, 500, rep)
                ).injection.n_attempts for rep in 1:10]
    @test all(>=(1), attempts)
    @test sum(attempts) > 10

    # An accepted result is reproducible on its own from its recorded seed.
    r9 = run_sel1_accepted(birth, death9,
             Sel1Params(1.0, 0.9, 5.0, 2.0, 1_000, :gamma, 0.1, 500, 1))
    again = run_sel1_once(birth, death9, 1_000, 0.1, 500,
                          MersenneTwister(r9.injection.seed))
    @test count(>(1.0), fitness_per_cell(again.pop)) == r9.injection.driver_clone_size

    # An unsatisfiable cell errors rather than returning junk: N_critic above
    # N_target means the hook can never fire, so nothing is ever accepted.
    p_bad = Sel1Params(1.0, 0.0, 5.0, 2.0, 500, :gamma, 0.5, 5_000, 1)
    @test_throws ErrorException run_sel1_accepted(birth, nodeath, p_bad;
                                                  max_attempts = 3)
end

@testset "clone size responds to s and to N_critic" begin
    birth = f -> Gamma(5.0, 1.0 / (5.0 * f))
    death = _ -> Gamma(5.0, 1.0 / (5.0 * 0.5))

    mean_clone(s, nc) = mean(
        run_sel1_accepted(birth, death,
            Sel1Params(1.0, 0.5, 5.0, 2.0, 1_000, :gamma, s, nc, rep)
        ).injection.driver_clone_size for rep in 1:5)

    # A stronger driver sweeps further by the time the population reaches N_target.
    @test mean_clone(2.0, 16) > mean_clone(0.1, 16)
    # An earlier driver has more generations to expand in.
    @test mean_clone(1.0, 2) > mean_clone(1.0, 500)

    # Retries get rarer as the driver gets stronger (establishment is easier).
    # Establishment rate is measured at d = 0.9, not at the d = 0.5 used above. At
    # d = 0.5 the driver establishes easily for every s (mean attempts 1.1-1.5, against
    # a floor of 1), so the comparison has almost no dynamic range and has no reliable
    # direction at N_critic = 500: it measured 1.22 vs 1.16 under one seed set and
    # 1.2 vs 1.5 under another, i.e. it changes sign with the seeds.
    # At d = 0.9 the same comparison spans 9.26 vs 5.64, a factor of 1.64. N_critic = 500
    # holds population-extinction risk constant between the two s values, so the
    # comparison isolates driver establishment. Measured 2026-09-02.
    death9 = _ -> Gamma(5.0, 1.0 / (5.0 * 0.9))
    mean_attempts(s) = mean(
        run_sel1_accepted(birth, death9,
            Sel1Params(1.0, 0.9, 5.0, 2.0, 1_000, :gamma, s, 500, rep)
        ).injection.n_attempts for rep in 1:50)
    @test mean_attempts(2.0) < mean_attempts(0.1)
end

end
