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
# Selection scenario 2 — deterministic (Dirac) division timing, the k → ∞ limit.
#
# Every mutation carries a fitness effect. Each daughter draws j ~ Poisson(ν)
# mutations at division; each one draws X ~ Gamma(a, 1/a) (mean 1) and updates
#
#     f ← min( f · (1 + s · X · (1 − f/M)),  M )
#
# See `growth_sel2_gamma.jl` for the reasoning behind ν = 1.0 (rather than the
# neutral runs' 2.0) and behind M = 10.0.
#
# Unlike the neutral deterministic runs, the division time is fitness-scaled: a cell
# of fitness f divides exactly 1/(b·f) after its birth. Fitness heterogeneity
# therefore breaks the lockstep, and the tree is no longer a perfect binary tree —
# leaf depths spread out where the neutral run has them all equal to log2(N). That
# loss of balance is the whole point of the panel: selection acting on tree shape
# with no timing noise at all. At s = 0 the lockstep and the perfect binary tree
# return exactly (verified in `verify_sel2.jl`).
#
# N_VALUES are powers of two, so the run stops at exactly N_target in both cases;
# the neutral runs' overshoot-to-the-next-power-of-two only shows up for a target
# that is not itself a power of two.
#
# d = 0 only, and N_VALUES are powers of two, both matching
# `analysis/sim_runs/neutral/growth_neutral_deterministic.jl` so each output file has
# an exactly paired neutral file. With no death there is no extinction, so
# `n_restarts` is always 0 here.

const N_SIMS       = 200
const b            = 1.0   # cell of fitness f divides at exactly 1/(b·f) after birth
const ν            = 1.0   # mutations per daughter per division (all drivers)
const M            = 10.0  # fitness cap
const EFFECT_SHAPE = 2.0   # a of Gamma(a, 1/a); CV = 1/√a ≈ 0.71
const TRAJ_DT      = 0.02
const N_VALUES     = [1_024, 16_384]
const S_VALUES     = [0.05, 0.10, 0.15, 0.20]

include(joinpath(HELPERS, "types_selection2.jl"))
include(joinpath(HELPERS, "selection2.jl"))

# ── Main loop ─────────────────────────────────────────────────────────────────

birth = f -> Dirac(1.0 / (b * f))
death = _ -> Dirac(1e8)   # no death; same sentinel as the neutral deterministic runs

for N_target in N_VALUES, s in S_VALUES
    println("Running $N_SIMS simulations  deterministic  N=$N_target  s=$s  M=$M ...")

    results = Vector{Sel2SimResult}(undef, N_SIMS)
    Threads.@threads for rep in 1:N_SIMS
        # gamma_shape = Inf and d = 0.0, matching the neutral convention.
        params = Sel2Params(b, 0.0, Inf, ν, N_target, :deterministic, s, M,
                            EFFECT_SHAPE, rep,
                            sel2_seed(:deterministic, N_target, 0.0, s, rep))
        results[rep] = run_sel2_once(birth, death, params; trajectory_dt = TRAJ_DT)
    end

    # No d or k in the filename, matching neutral_deterministic_N1024.jls.
    outfile = joinpath(DATA, "raw", "selection_2",
                       "deterministic",
                       "sel2_deterministic_N$(N_target)_s$(s)_M$(M).jls")
    mkpath(dirname(outfile))
    serialize(outfile, results)
    println("  → $(outfile)")
end

println("All done.")
