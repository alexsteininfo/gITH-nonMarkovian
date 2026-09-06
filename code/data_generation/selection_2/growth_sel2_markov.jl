using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using Distributions
using Serialization
using Random

# ── Parameters ────────────────────────────────────────────────────────────────
# Selection scenario 2 — Markovian (exponential) division timing, i.e. the k = 1
# Gamma / standard Gillespie birth-death model. This is the reference against which
# the Gamma and deterministic timing models are read.
#
# Every mutation carries a fitness effect. Each daughter draws j ~ Poisson(ν)
# mutations at division; each one draws X ~ Gamma(a, 1/a) (mean 1) and updates
#
#     f ← min( f · (1 + s · X · (1 − f/M)),  M )
#
# See `growth_sel2_gamma.jl` for the reasoning behind ν = 1.0 (rather than the
# neutral runs' 2.0) and behind M = 10.0.
#
# Everything else — b, the N/d grid, the fitness-independent death — is held identical
# to `analysis/sim_runs/neutral/growth_neutral_markov.jl`, so each output file has an
# exactly paired neutral file.

const N_SIMS       = 200
const b            = 1.0   # baseline birth rate at f = 1; mean division time 1/(b·f)
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
    println("Running $N_SIMS simulations  markov  N=$N_target  d=$d  s=$s  M=$M ...")

    # Birth is fitness-scaled (mean 1/(b·f)); death is fitness-independent.
    # d = 0 uses the neutral runs' 1e8 scale sentinel so death never fires.
    death_scale = d > 0.0 ? 1.0 / d : 1e8
    birth = f -> Exponential(1.0 / (b * f))
    death = _ -> Exponential(death_scale)

    results = Vector{Sel2SimResult}(undef, N_SIMS)
    Threads.@threads for rep in 1:N_SIMS
        # gamma_shape = 1.0 for the markov model, matching the neutral convention.
        params = Sel2Params(b, d, 1.0, ν, N_target, :markov, s, M, EFFECT_SHAPE, rep,
                            sel2_seed(:markov, N_target, d, s, rep))
        results[rep] = run_sel2_once(birth, death, params; trajectory_dt = TRAJ_DT)
    end

    # No k in the filename, matching neutral_markov_N10000_d0.5.jls.
    outfile = joinpath(@__DIR__, "..", "..", "..", "data", "raw", "selection_2", "markov",
                       "sel2_markov_N$(N_target)_d$(d)_s$(s)_M$(M).jls")
    mkpath(dirname(outfile))
    serialize(outfile, results)
    println("  → $(outfile)")
end

println("All done.")
