using Pkg
const ROOT = dirname(dirname(dirname(dirname(@__DIR__))))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using Distributions
using Serialization
using Statistics
using Random

# ── Parameters ────────────────────────────────────────────────────────────────
# Dirac waiting times: a cell of fitness f divides at exactly 1/(b·f) after birth.
# No death (d = 0), so no run is ever lost to extinction and the driver clone can
# never die out — every cell is accepted on its first attempt.
# Background process identical to growth_neutral_deterministic.jl.
#
# All cells divide in lockstep, so N_critic matters only through its generation
# ceil(log2(N_critic)); the ratio-2 grid gives one N_critic per generation, and
# N_target = 1024 lands exactly on 2^0 … 2^9.

const b        = 1.0
const k        = Inf   # recorded as gamma_shape: the k → ∞ limit
const nu       = 2.0
const N_VALUES = [1_024, 16_384]
const S_VALUES = collect(0.1:0.1:2.0)
const N_REPS   = 5

include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "selection1.jl"))

# ── Main loop ─────────────────────────────────────────────────────────────────

for N_target in N_VALUES
    birth_dist = f -> Dirac(1.0 / (b * f))
    death_dist = _ -> Dirac(1e8)
    N_CRITIC   = ncritic_grid(N_target)

    for s in S_VALUES
        # The sweep is resumable: a shard is skipped once its output file exists, so a
        # crash mid-sweep restarts only from the first missing shard. Delete a shard's
        # .jls file to force it to be regenerated.
        outfile = joinpath(DATA, "raw", "selection_1",
                           "deterministic", "sel1_deterministic_N$(N_target)_s$(s).jls")
        if isfile(outfile)
            println("  → $(outfile) exists, skipping")
            continue
        end

        println("deterministic  N=$N_target  s=$s ...")

        cells   = [(nc, rep) for nc in N_CRITIC for rep in 1:N_REPS]
        results = Vector{Sel1SimResult}(undef, length(cells))
        Threads.@threads for i in eachindex(cells)
            nc, rep = cells[i]
            results[i] = run_sel1_accepted(birth_dist, death_dist,
                Sel1Params(b, 0.0, k, nu, N_target, :deterministic, s, nc, rep))
        end

        mkpath(dirname(outfile))
        serialize(outfile, results)
        println("  → $(outfile)")
    end
end

println("All done.")
