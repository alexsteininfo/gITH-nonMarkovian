using Pkg
const ROOT = dirname(dirname(dirname(dirname(@__DIR__))))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using AbstractTrees
using Serialization

include(joinpath(HELPERS, "types_selection1.jl"))

# ── Paths ─────────────────────────────────────────────────────────────────────

const RAW  = joinpath(DATA, "raw",       "selection_1")
const PROC = joinpath(DATA, "processed", "selection_1")

# ── Processing ────────────────────────────────────────────────────────────────

# Stage 2 for selection scenario 1: one driver of strength s injected at N_critic.
# Same four quantities as `process_neutral.jl`, plus the injection record.
#
#   params       — Sel1Params (carries s, N_critic, rep; nu = 2.0 as in neutral)
#   mut_per_cell — per-cell neutral mutation burden (one count per alive leaf)
#   sfs          — site frequency spectrum, sfs[k] = # mutations carried by k cells
#   leaf_depths  — divisions from root to each leaf
#   injection    — Sel1Injection: t_inject, N_at_inject, driver_cell_id,
#                  driver_clone_size, n_attempts, seed
#
# `injection` is saved so the plotting stage never has to reopen the 6.9 GB of raw
# trees for what are per-simulation scalars — `driver_clone_size` is the headline
# observable of this scenario and `1/mean(n_attempts)` its establishment probability.
#
# The driver does *not* increment `mutations` (see `helpers/selection1.jl`), so the
# neutral ν = 2.0 channel is untouched and `sfs`/`mut_per_cell` stay directly
# comparable to the paired `neutral/` file. The clone is identified by fitness > 1,
# recorded as a count in `injection.driver_clone_size`.
#
# Each shard holds 50 sims: 10 log-spaced N_critic × 5 reps, all at one (model, N, d, s).

const SUBDIRS = ("params", "mut_per_cell", "sfs", "leaf_depths", "injection")

function process_file(raw_path::String, proc_dir::String, label::String)
    stem = splitext(basename(raw_path))[1]

    # Resumable, like the sweep itself: delete the processed file to force a redo.
    if all(isfile(joinpath(proc_dir, sub, stem * ".jls")) for sub in SUBDIRS)
        println("  Skipping (already processed): $(basename(raw_path))")
        return
    end

    println("  Processing: $(basename(raw_path))")

    sims = deserialize(raw_path)::Vector{Sel1SimResult}

    all_params       = Sel1Params[]
    all_mut_per_cell = Vector{Int}[]
    all_sfs          = Vector{Int}[]
    all_leaf_depths  = Vector{Int}[]
    all_injection    = Sel1Injection[]

    # `run_sel1_accepted` only returns runs with a unique root, so nothing trees should
    # not occur here — kept for symmetry with the neutral and scenario-2 scripts, and
    # because a warning is cheaper than a silent misalignment if that ever changes.
    skipped = 0
    for sim in sims
        if isnothing(sim.tree_root)
            skipped += 1
            continue
        end
        root = sim.tree_root

        # Actual leaf count: deterministic timing overshoots N_target once the driver
        # divides off-phase, and every model stops at the first popsize >= N_target.
        actual_N = length(collect(Leaves(root)))

        push!(all_params,       sim.params)
        push!(all_mut_per_cell, mutations_per_cell(root))
        push!(all_sfs,          sitefrequencyspectrum(root, actual_N))
        push!(all_leaf_depths,  leaf_depths(root))
        push!(all_injection,    sim.injection)
    end

    skipped > 0 && @warn "  $label: skipped $skipped / $(length(sims)) sims with nothing tree"

    for (subdir, data) in (
        ("params",        all_params),
        ("mut_per_cell",  all_mut_per_cell),
        ("sfs",           all_sfs),
        ("leaf_depths",   all_leaf_depths),
        ("injection",     all_injection),
    )
        outdir  = joinpath(proc_dir, subdir)
        mkpath(outdir)
        serialize(joinpath(outdir, stem * ".jls"), data)
    end

    println("    → saved $(length(all_params)) results to $proc_dir")
end

# The grids below mirror `analysis/sim_runs/selection_1/growth_sel1_*.jl` exactly,
# including `collect(0.1:0.1:2.0)` rather than a literal list, so the interpolated
# filenames match the shards on disk character for character.

const S_VALUES = collect(0.1:0.1:2.0)

# ── Gamma shards (120): N ∈ {1000, 10000} × d ∈ {0.0, 0.5, 0.9} × 20 s, k = 5.0 ──

println("\n── Gamma model ────────────────────────────────────────────────────────")
gamma_raw  = joinpath(RAW,  "gamma")
gamma_proc = joinpath(PROC, "gamma")

for N in (1_000, 10_000), d in (0.0, 0.5, 0.9), s in S_VALUES
    fname = "sel1_gamma_N$(N)_d$(d)_k5.0_s$(s).jls"
    raw_path = joinpath(gamma_raw, fname)
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, gamma_proc, "gamma N=$N d=$d s=$s")
end

# ── Markov shards (120): N ∈ {1000, 10000} × d ∈ {0.0, 0.5, 0.9} × 20 s ─────────

println("\n── Markov model ───────────────────────────────────────────────────────")
markov_raw  = joinpath(RAW,  "markov")
markov_proc = joinpath(PROC, "markov")

for N in (1_000, 10_000), d in (0.0, 0.5, 0.9), s in S_VALUES
    fname = "sel1_markov_N$(N)_d$(d)_s$(s).jls"
    raw_path = joinpath(markov_raw, fname)
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, markov_proc, "markov N=$N d=$d s=$s")
end

# ── Deterministic shards (40): N ∈ {1024, 16384} × 20 s, d = 0 ─────────────────

println("\n── Deterministic model ────────────────────────────────────────────────")
det_raw  = joinpath(RAW,  "deterministic")
det_proc = joinpath(PROC, "deterministic")

for N in (1_024, 16_384), s in S_VALUES
    fname = "sel1_deterministic_N$(N)_s$(s).jls"
    raw_path = joinpath(det_raw, fname)
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, det_proc, "deterministic N=$N s=$s")
end

println("\nAll done.")
