using Pkg
Pkg.activate(dirname(dirname(@__DIR__)))

using MutationLoadDynamics
using AbstractTrees
using Serialization

include(joinpath(dirname(@__DIR__), "helpers", "types_selection2.jl"))

# ── Paths ─────────────────────────────────────────────────────────────────────

const RAW  = joinpath(@__DIR__, "..", "..", "data", "raw",       "selection_2")
const PROC = joinpath(@__DIR__, "..", "..", "data", "processed", "selection_2")

# ── Processing ────────────────────────────────────────────────────────────────

# Stage 2 for selection scenario 2: every mutation is a driver, capped at M = 10.
# Same four quantities as `process_neutral.jl`, plus fitness and the restart count.
#
#   params       — Sel2Params (carries s, M, effect_shape, rep, seed)
#   mut_per_cell — per-cell mutation burden (one count per alive leaf)
#   sfs          — site frequency spectrum, sfs[k] = # mutations carried by k cells
#   leaf_depths  — divisions from root to each leaf
#   leaf_fitness — final fitness of each alive leaf
#   n_restarts   — extinctions before the attempt that reached N_target
#
# Two differences from the neutral runs that matter when reading these arrays:
#
#   ν = 1.0, not 2.0. SFS *shapes* still overlay the neutral runs after a factor-ν
#   rescale; the single-cell burden distribution does not, so an scMB comparison
#   against `neutral/` needs neutral runs at ν = 1.
#
#   There is no separate neutral channel — the mutations counted in `sfs` and
#   `mut_per_cell` are the same ones carrying the fitness effects. `mut_per_cell` is
#   therefore also the driver count per cell, and it is co-indexed with
#   `leaf_fitness` (both iterate the alive leaves in `Leaves` order), so burden and
#   fitness can be paired cell by cell. `leaf_depths` uses a different traversal and
#   is a pooled distribution only.
#
# `leaf_fitness` and `n_restarts` are saved so the plotting stage never has to reopen
# the 4.9 GB of raw trees: fitness is this scenario's central observable, and
# `1/(1 + mean(n_restarts))` estimates founder survival probability.
#
# Each file holds 200 sims at one (model, N_target, d, s).

const SUBDIRS = ("params", "mut_per_cell", "sfs", "leaf_depths", "leaf_fitness", "n_restarts")

function process_file(raw_path::String, proc_dir::String, label::String)
    stem = splitext(basename(raw_path))[1]

    # Resumable: delete the processed file to force a redo.
    if all(isfile(joinpath(proc_dir, sub, stem * ".jls")) for sub in SUBDIRS)
        println("  Skipping (already processed): $(basename(raw_path))")
        return
    end

    println("  Processing: $(basename(raw_path))")

    sims = deserialize(raw_path)::Vector{Sel2SimResult}

    all_params       = Sel2Params[]
    all_mut_per_cell = Vector{Int}[]
    all_sfs          = Vector{Int}[]
    all_leaf_depths  = Vector{Int}[]
    all_leaf_fitness = Vector{Float64}[]
    all_n_restarts   = Int[]

    # `run_sel2_once` documents that `getsingleroot` may find no unique root; skip
    # those, exactly as the neutral script does.
    skipped = 0
    for sim in sims
        if isnothing(sim.tree_root)
            skipped += 1
            continue
        end
        root = sim.tree_root

        # Actual leaf count: the run stops at the first popsize >= N_target, and
        # deterministic timing overshoots once fitness differences desynchronise it.
        actual_N = length(collect(Leaves(root)))

        push!(all_params,       sim.params)
        push!(all_mut_per_cell, mutations_per_cell(root))
        push!(all_sfs,          sitefrequencyspectrum(root, actual_N))
        push!(all_leaf_depths,  leaf_depths(root))
        push!(all_leaf_fitness, leaf_fitness(root))
        push!(all_n_restarts,   sim.n_restarts)
    end

    skipped > 0 && @warn "  $label: skipped $skipped / $(length(sims)) sims with nothing tree"

    for (subdir, data) in (
        ("params",        all_params),
        ("mut_per_cell",  all_mut_per_cell),
        ("sfs",           all_sfs),
        ("leaf_depths",   all_leaf_depths),
        ("leaf_fitness",  all_leaf_fitness),
        ("n_restarts",    all_n_restarts),
    )
        outdir  = joinpath(proc_dir, subdir)
        mkpath(outdir)
        serialize(joinpath(outdir, stem * ".jls"), data)
    end

    println("    → saved $(length(all_params)) results to $proc_dir")
end

# The grids below mirror `analysis/sim_runs/selection_2/growth_sel2_*.jl` exactly, so
# the interpolated filenames match the files on disk character for character —
# note `0.10` and `0.20` interpolate as `s0.1` and `s0.2`.

const S_VALUES = [0.05, 0.10, 0.15, 0.20]
const M_CAP    = 10.0

# ── Gamma files (24): N ∈ {1000, 10000} × d ∈ {0.0, 0.5, 0.9} × 4 s, k = 5.0 ────

println("\n── Gamma model ────────────────────────────────────────────────────────")
gamma_raw  = joinpath(RAW,  "gamma")
gamma_proc = joinpath(PROC, "gamma")

for N in (1_000, 10_000), d in (0.0, 0.5, 0.9), s in S_VALUES
    fname = "sel2_gamma_N$(N)_d$(d)_k5.0_s$(s)_M$(M_CAP).jls"
    raw_path = joinpath(gamma_raw, fname)
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, gamma_proc, "gamma N=$N d=$d s=$s")
end

# ── Markov files (24): N ∈ {1000, 10000} × d ∈ {0.0, 0.5, 0.9} × 4 s ───────────

println("\n── Markov model ───────────────────────────────────────────────────────")
markov_raw  = joinpath(RAW,  "markov")
markov_proc = joinpath(PROC, "markov")

for N in (1_000, 10_000), d in (0.0, 0.5, 0.9), s in S_VALUES
    fname = "sel2_markov_N$(N)_d$(d)_s$(s)_M$(M_CAP).jls"
    raw_path = joinpath(markov_raw, fname)
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, markov_proc, "markov N=$N d=$d s=$s")
end

# ── Deterministic files (8): N ∈ {1024, 16384} × 4 s, d = 0 ────────────────────

println("\n── Deterministic model ────────────────────────────────────────────────")
det_raw  = joinpath(RAW,  "deterministic")
det_proc = joinpath(PROC, "deterministic")

for N in (1_024, 16_384), s in S_VALUES
    fname = "sel2_deterministic_N$(N)_s$(s)_M$(M_CAP).jls"
    raw_path = joinpath(det_raw, fname)
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, det_proc, "deterministic N=$N s=$s")
end

println("\nAll done.")
