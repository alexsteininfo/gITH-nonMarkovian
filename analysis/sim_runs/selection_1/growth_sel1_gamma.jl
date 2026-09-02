using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Distributions
using Serialization
using Statistics
using Random

# ── Parameters ────────────────────────────────────────────────────────────────
# Background process identical to analysis/sim_runs/neutral/growth_neutral_gamma.jl.
# The only difference is the single injected driver; see
# docs/superpowers/specs/2026-09-01-selection-scenario-1-design.md.

const b        = 1.0   # birth rate; mean division time = 1/(b·f)
const k        = 5.0   # Gamma shape; CV = 1/√k ≈ 0.45 (cancer-cell regime)
const nu       = 2.0   # neutral mutations per daughter per division
const N_VALUES = [1_000, 10_000]
const D_VALUES = [0.0, 0.5, 0.9]
const S_VALUES = collect(0.1:0.1:2.0)
const N_REPS   = 5

const HELPERS = joinpath(dirname(dirname(@__DIR__)), "helpers")
include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "selection1.jl"))

# ── Main loop ─────────────────────────────────────────────────────────────────
# One output file per s (50 sims: 10 N_critic × 5 reps). Sharding by s keeps the
# worst file near 95 MB instead of 1.9 GB for a whole parameter set.

for N_target in N_VALUES, d in D_VALUES
    death_scale = d > 0.0 ? 1.0 / (k * d) : 1e8
    birth_dist  = f -> Gamma(k, 1.0 / (k * b * f))
    death_dist  = _ -> Gamma(k, death_scale)
    N_CRITIC    = ncritic_grid(N_target)

    for s in S_VALUES
        println("gamma  N=$N_target  d=$d  s=$s  ($(length(N_CRITIC)) N_critic × $N_REPS reps) ...")

        cells   = [(nc, rep) for nc in N_CRITIC for rep in 1:N_REPS]
        results = Vector{Sel1SimResult}(undef, length(cells))
        Threads.@threads for i in eachindex(cells)
            nc, rep = cells[i]
            results[i] = run_sel1_accepted(birth_dist, death_dist,
                Sel1Params(b, d, k, nu, N_target, :gamma, s, nc, rep))
        end

        outfile = joinpath(@__DIR__, "..", "..", "..", "data", "raw", "selection_1",
                           "gamma", "sel1_gamma_N$(N_target)_d$(d)_k$(k)_s$(s).jls")
        mkpath(dirname(outfile))
        serialize(outfile, results)
        println("  → $(outfile)  (mean attempts $(round(mean(r.injection.n_attempts for r in results), digits = 2)))")
    end
end

println("All done.")
