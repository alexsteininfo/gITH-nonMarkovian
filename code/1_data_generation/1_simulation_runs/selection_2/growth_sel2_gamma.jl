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
# Selection scenario 2 — Gamma (non-Markovian) division timing, k = 5.
#
# Every mutation carries a fitness effect. Each daughter draws j ~ Poisson(ν)
# mutations at division; each one draws X ~ Gamma(a, 1/a) (mean 1) and updates
#
#     f ← min( f · (1 + s · X · (1 − f/M)),  M )
#
# so `s` is the mean per-mutation selection coefficient and `a = EFFECT_SHAPE` sets
# only its spread (CV = 1/√a). There is no separate neutral channel: the same
# mutations that carry the fitness effects are the ones `sitefrequencyspectrum` and
# `mutations_per_cell` count.
#
# ν = 1.0, not the neutral runs' 2.0. At ν = 2 a lineage accumulates ~27 mutations by
# N = 10000 and the population piles onto the cap, collapsing fitness *differences*
# into an effectively neutral process with a faster clock. At ν = 1 the worst corner
# (s = 0.2, N = 10000, ~14 mutations per lineage) reaches f ≈ 6.1 of M = 10 with a
# 5-95 percentile spread of [4.2, 8.1], so all four s values stay graded.
# Consequence for plotting: SFS shapes still overlay the ν = 2 neutral runs after a
# factor-ν rescale, but the single-cell burden distribution does not — a fair scMB
# overlay needs neutral runs at ν = 1.
#
# Everything else — b, k, the N/d grid, the fitness-independent death — is held
# identical to `analysis/sim_runs/neutral/growth_neutral_gamma.jl`, so each output
# file has an exactly paired neutral file.

const N_SIMS       = 200
const b            = 1.0   # baseline birth rate at f = 1; mean division time 1/(b·f)
const k            = 5.0   # timing Gamma shape; CV = 1/√k ≈ 0.45 (cancer-cell regime)
const ν            = 1.0   # mutations per daughter per division (all drivers)
const M            = 10.0  # fitness cap
const EFFECT_SHAPE = 2.0   # a of Gamma(a, 1/a); CV = 1/√a ≈ 0.71
const TRAJ_DT      = 0.02
const N_VALUES     = [1_000, 10_000]
const D_VALUES     = [0.0, 0.5, 0.9]
const S_VALUES     = [0.05, 0.10, 0.15, 0.20]

include(joinpath(HELPERS, "types_selection2.jl"))
include(joinpath(HELPERS, "selection2.jl"))

# ── Main loop ─────────────────────────────────────────────────────────────────

for N_target in N_VALUES, d in D_VALUES, s in S_VALUES
    println("Running $N_SIMS simulations  gamma  N=$N_target  d=$d  k=$k  s=$s  M=$M ...")

    # Birth is fitness-scaled (mean 1/(b·f)); death is fitness-independent.
    # d = 0 uses the neutral runs' 1e8 scale sentinel so death never fires.
    death_scale = d > 0.0 ? 1.0 / (k * d) : 1e8
    birth = f -> Gamma(k, 1.0 / (k * b * f))
    death = _ -> Gamma(k, death_scale)

    results = Vector{Sel2SimResult}(undef, N_SIMS)
    Threads.@threads for rep in 1:N_SIMS
        params = Sel2Params(b, d, k, ν, N_target, :gamma, s, M, EFFECT_SHAPE, rep,
                            sel2_seed(:gamma, N_target, d, s, rep))
        results[rep] = run_sel2_once(birth, death, params; trajectory_dt = TRAJ_DT)
    end

    outfile = joinpath(DATA, "raw", "selection_2", "gamma",
                       "sel2_gamma_N$(N_target)_d$(d)_k$(k)_s$(s)_M$(M).jls")
    mkpath(dirname(outfile))
    serialize(outfile, results)
    println("  → $(outfile)")
end

println("All done.")
