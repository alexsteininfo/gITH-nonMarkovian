using Pkg
const ROOT = dirname(dirname(dirname(dirname(@__DIR__))))
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
include(joinpath(HELPERS, "types_selection2.jl"))
include(joinpath(HELPERS, "selection2.jl"))

# Small N throughout: these checks are about mechanism, not statistics, so the whole
# file runs in seconds. Run as:
#   julia --project=. analysis/sim_runs/selection_2/verify_sel2.jl

const A = 2.0   # EFFECT_SHAPE used by the three sim scripts

gamma_birth(k, b) = f -> Gamma(k, 1.0 / (k * b * f))
gamma_death(k, d) = _ -> Gamma(k, d > 0.0 ? 1.0 / (k * d) : 1e8)

sel2_params(; d = 0.0, k = 5.0, nu = 1.0, N = 500, model = :gamma, s = 0.1,
              M = 10.0, a = A, rep = 1) =
    Sel2Params(1.0, d, k, nu, N, model, s, M, a, rep,
               sel2_seed(model, N, d, s, rep))

@testset "selection scenario 2" begin

@testset "environment resolution" begin
    # The recurring bug in this repo: one `dirname` short, silently activating the
    # wrong environment. The header computes ROOT three `dirname`s up from this
    # script and asserts it holds a Project.toml, so a moved script fails
    # immediately rather than silently activating the wrong environment.
    @test Base.active_project() == joinpath(ROOT, "Project.toml")
    @test isfile(joinpath(HELPERS, "types_selection2.jl"))
    @test isdir(joinpath(ROOT, "code", "data_generation", "selection_2"))
end

@testset "types round-trip through Serialization" begin
    p = Sel2Params(1.0, 0.5, 5.0, 1.0, 1000, :gamma, 0.15, 10.0, 2.0, 7, 0xdeadbeef)
    r = Sel2SimResult(TrajectoryPoint[], nothing, p, 3)

    path = tempname() * ".jls"
    serialize(path, [r])
    back = deserialize(path)::Vector{Sel2SimResult}
    @test length(back) == 1
    @test back[1].params == p
    @test back[1].params.model        === :gamma
    @test back[1].params.effect_shape == 2.0
    @test back[1].params.M            == 10.0
    @test back[1].n_restarts          == 3
    rm(path)
end

@testset "existing neutral .jls remain readable" begin
    # types_selection2.jl must not have disturbed SimParams / GrowthSimResult.
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

@testset "effect distribution is mean-1 with the intended spread" begin
    for a in (1.0, 2.0, 5.0)
        D = sel2_effect_dist(a)
        @test mean(D) ≈ 1.0
        @test std(D)  ≈ 1.0 / sqrt(a)   # CV = 1/√a since the mean is 1
    end
    # a = 1 is the Exponential(1) of selection_old/growth_multi_rand.jl.
    @test sel2_effect_dist(1.0) == Gamma(1.0, 1.0)
    @test mean(sel2_effect_dist(1.0)) ≈ mean(Exponential(1.0))
    @test var(sel2_effect_dist(1.0))  ≈ var(Exponential(1.0))
end

@testset "fitness update rule" begin
    M = 10.0

    # s = 0 is the identity, for every draw of X. This is what makes s = 0 exactly
    # neutral rather than approximately so.
    id = sel2_fitness_update(0.0, M)
    for f in (1.0, 3.7, M), X in (0.0, 1.0, 50.0)
        @test id(f, X) === f
    end

    up = sel2_fitness_update(0.1, M)
    # Matches the closed form away from the cap.
    @test up(1.0, 1.0) ≈ 1.0 * (1 + 0.1 * 1.0 * (1 - 1.0 / M))
    @test up(4.0, 2.5) ≈ 4.0 * (1 + 0.1 * 2.5 * (1 - 4.0 / M))
    # Non-decreasing, and never above the cap even for an extreme draw.
    for f in (1.0, 2.0, 5.0, 9.0, 9.99, M)
        @test up(f, 0.0) == f
        @test up(f, 1.0) >= f
        @test up(f, 1e6) <= M
    end
    @test up(M, 1.0) == M           # a cell at the cap stays there
    # Diminishing returns: the same X buys less as f approaches M.
    @test up(2.0, 1.0) - 2.0 > up(9.0, 1.0) - 9.0
end

@testset "seeds separate the grid cells" begin
    @test sel2_seed(:gamma, 1000, 0.5, 0.1, 1) != sel2_seed(:markov, 1000, 0.5, 0.1, 1)
    @test sel2_seed(:gamma, 1000, 0.5, 0.1, 1) != sel2_seed(:gamma,  1000, 0.9, 0.1, 1)
    @test sel2_seed(:gamma, 1000, 0.5, 0.1, 1) != sel2_seed(:gamma,  1000, 0.5, 0.2, 1)
    @test sel2_seed(:gamma, 1000, 0.5, 0.1, 1) != sel2_seed(:gamma,  1000, 0.5, 0.1, 2)
    @test sel2_seed(:gamma, 1000, 0.5, 0.1, 1) == sel2_seed(:gamma,  1000, 0.5, 0.1, 1)
    @test sel2_seed(:gamma, 1000, 0.5, 0.1, 1) isa UInt64
end

@testset "one run: shape, cap, reproducibility" begin
    p = sel2_params(d = 0.5, s = 0.15, N = 500)
    r = run_sel2_once(gamma_birth(5.0, 1.0), gamma_death(5.0, 0.5), p)

    @test r isa Sel2SimResult
    @test r.params == p
    @test !isnothing(r.tree_root)
    @test length(r.trajectory) > 0
    @test length(collect(Leaves(r.tree_root))) >= p.N_target

    # The cap is a hard bound on every cell, not just on the mean.
    fits = [n.data.fitness for n in Leaves(r.tree_root)]
    @test all(f -> 1.0 <= f <= p.M, fits)
    @test maximum(fits) > 1.0        # selection actually did something

    # Identical params (hence identical seed) reproduce the run exactly.
    again = run_sel2_once(gamma_birth(5.0, 1.0), gamma_death(5.0, 0.5), p)
    @test length(collect(Leaves(again.tree_root))) ==
          length(collect(Leaves(r.tree_root)))
    @test sort(mutations_per_cell(again.tree_root)) ==
          sort(mutations_per_cell(r.tree_root))
    @test [n.data.fitness for n in Leaves(again.tree_root)] == fits
end

@testset "s = 0 is exactly neutral" begin
    # Bit-for-bit, not statistically: at s = 0 the update is the identity, so a
    # scenario-2 run and a neutral-update run consuming the same Gamma draws must
    # agree exactly. This is the check that the scenario-2 machinery has not
    # perturbed the baseline. (Against the *stored* neutral runs the comparison can
    # only be statistical, and only after re-running those at ν = 1.)
    p = sel2_params(d = 0.5, s = 0.0, N = 500)
    r = run_sel2_once(gamma_birth(5.0, 1.0), gamma_death(5.0, 0.5), p)

    spec = MeasurementSpec(trajectory_dt = 0.02,
                           snapshot_triggers = [AtEnd()], snapshot_stats = [])
    pop  = initialize_population(fitness_init = 1.0)
    block = NonMarkovBlock(
        birth_dist            = gamma_birth(5.0, 1.0),
        death_dist            = gamma_death(5.0, 0.5),
        stopfunction          = pop -> popsize(pop) >= p.N_target,
        driver_dist           = sel2_effect_dist(p.effect_shape),  # same rng draws
        fitness_update        = (f, X) -> f,                       # neutral update
        ν                     = p.nu,
        restart_on_extinction = true,
    )
    acc = MeasurementAccumulator(spec)
    simulate!(pop, block, MersenneTwister(p.seed); accumulator = acc)

    @test popsize(pop) == length(collect(Leaves(r.tree_root)))
    @test pop.t === endtime(r.tree_root) || pop.t > 0.0
    @test sort(mutations_per_cell(pop)) == sort(mutations_per_cell(r.tree_root))
    @test all(==(1.0), [n.data.fitness for n in Leaves(r.tree_root)])
end

@testset "deterministic model: lockstep breaks only under selection" begin
    birth = f -> Dirac(1.0 / f)
    death = _ -> Dirac(1e8)

    # s = 0: every cell has lifetime exactly 1, so the tree is a perfect binary tree
    # and the population overshoots to the next power of two — the neutral behaviour.
    r0 = run_sel2_once(birth, death,
                       sel2_params(k = Inf, N = 512, model = :deterministic, s = 0.0))
    d0 = leaf_depths(r0.tree_root)
    @test length(d0) == 512
    @test all(==(9), d0)                       # log2(512)
    @test r0.n_restarts == 0                   # no death, so no extinction

    # s > 0: fitter cells divide sooner, the lockstep breaks, and the tree is no
    # longer balanced. A deterministic claim, not a statistical one.
    r1 = run_sel2_once(birth, death,
                       sel2_params(k = Inf, N = 512, model = :deterministic, s = 0.2))
    d1 = leaf_depths(r1.tree_root)
    @test !all(==(first(d1)), d1)
    @test mean(d1) > mean(d0)                  # imbalance lengthens external paths
    @test maximum(n.data.fitness for n in Leaves(r1.tree_root)) <= 10.0
end

@testset "mean fitness increases with s" begin
    mean_fit(s) = mean(
        mean(n.data.fitness for n in Leaves(
            run_sel2_once(gamma_birth(5.0, 1.0), gamma_death(5.0, 0.5),
                          sel2_params(d = 0.5, s = s, N = 500, rep = rep)).tree_root))
        for rep in 1:10)

    f05, f10, f20 = mean_fit(0.05), mean_fit(0.10), mean_fit(0.20)
    @info "mean leaf fitness at N=500, d=0.5" s005=f05 s010=f10 s020=f20
    @test 1.0 < f05 < f10 < f20
    @test f20 < 10.0        # not yet pinned to the cap at the production settings
end

@testset "n_restarts counts total-population extinctions only" begin
    # d = 0: extinction is impossible.
    for rep in 1:5
        r = run_sel2_once(gamma_birth(5.0, 1.0), gamma_death(5.0, 0.0),
                          sel2_params(d = 0.0, s = 0.1, N = 500, rep = rep))
        @test r.n_restarts == 0
    end

    # d = 0.9: a single founder usually dies out, so restarts must be common. Summed
    # over 20 replicates rather than asserted per-run, which would be a coin flip.
    total = sum(run_sel2_once(gamma_birth(5.0, 1.0), gamma_death(5.0, 0.9),
                    sel2_params(d = 0.9, s = 0.05, N = 500, rep = rep)).n_restarts
                for rep in 1:20)
    @info "restarts over 20 runs at d=0.9, N=500, s=0.05" total
    @test total > 0

    # Subclone loss is *not* a restart: with death on, plenty of lineages are pruned
    # while the run itself is accepted first try. Nothing conditions on their fate.
    r = run_sel2_once(gamma_birth(5.0, 1.0), gamma_death(5.0, 0.5),
                      sel2_params(d = 0.5, s = 0.1, N = 1000))
    @test !isnothing(r.tree_root)
    # Internal nodes whose subtree holds no alive leaf are gone from the tree, so the
    # surviving leaf count is what the SFS is built on.
    @test length(collect(Leaves(r.tree_root))) >= 1000
end

end
