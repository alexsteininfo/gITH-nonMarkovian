using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using AbstractTrees
using Serialization

include(joinpath(HELPERS, "types.jl"))

# ── Paths ─────────────────────────────────────────────────────────────────────

const RAW  = joinpath(DATA, "raw",       "neutral")
const PROC = joinpath(DATA, "processed", "neutral")

# ── Processing ────────────────────────────────────────────────────────────────

# For each simulation result, extract:
#   params       — SimParams struct (nu=2.0, s=0.0)
#   mut_per_cell — actual per-cell mutation burden (Vector{Int}, one count per leaf cell)
#   sfs          — site frequency spectrum (Vector{Int64}, sfs[k] = # mutations in k cells)
#
# With ν=2.0, each division sprinkles Poisson(2.0) neutral mutations onto the
# branch; sitefrequencyspectrum and mutations_per_cell read node.data.mutations directly.

function process_file(raw_path::String, proc_dir::String, label::String)
    println("  Processing: $(basename(raw_path))")

    sims = deserialize(raw_path)::Vector{GrowthSimResult}

    # Collect per-simulation outputs, skipping any run whose tree is nothing.
    # With ν=2.0, node.data.mutations holds actual simulated neutral mutations;
    # mutations_per_cell and sitefrequencyspectrum read these directly.
    # leaf_depths are also saved for empirical Theorem 3 (avoids Poisson D_L≈1 approximation).
    all_params       = SimParams[]
    all_mut_per_cell = Vector{Int}[]
    all_sfs          = Vector{Int}[]
    all_leaf_depths  = Vector{Int}[]

    skipped = 0
    for sim in sims
        if isnothing(sim.tree_root)
            skipped += 1
            continue
        end
        root = sim.tree_root

        # Actual leaf count (deterministic model overshoots N_target to next power of 2).
        actual_N = length(collect(Leaves(root)))

        push!(all_params,       sim.params)
        push!(all_mut_per_cell, mutations_per_cell(root))
        push!(all_sfs,          sitefrequencyspectrum(root, actual_N))
        push!(all_leaf_depths,  leaf_depths(root))
    end

    skipped > 0 && @warn "  $label: skipped $skipped / $(length(sims)) sims with nothing tree"

    stem = splitext(basename(raw_path))[1]

    for (subdir, data) in (
        ("params",        all_params),
        ("mut_per_cell",  all_mut_per_cell),
        ("sfs",           all_sfs),
        ("leaf_depths",   all_leaf_depths),
    )
        outdir  = joinpath(proc_dir, subdir)
        mkpath(outdir)
        outfile = joinpath(outdir, stem * ".jls")
        serialize(outfile, data)
    end

    println("    → saved $(length(all_params)) results to $proc_dir")
end

# ── Gamma files (6): N ∈ {1000, 10000}, d ∈ {0.0, 0.5, 0.9}, k = 5.0 ────────

println("\n── Gamma model ────────────────────────────────────────────────────────")
gamma_raw  = joinpath(RAW,  "gamma")
gamma_proc = joinpath(PROC, "gamma")

for N in (1_000, 10_000), d in (0.0, 0.5, 0.9)
    fname = "neutral_gamma_N$(N)_d$(d)_k5.0.jls"
    raw_path = joinpath(gamma_raw, fname)
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, gamma_proc, "gamma N=$N d=$d")
end

# ── Markov files (6): N ∈ {1000, 10000}, d ∈ {0.0, 0.5, 0.9} ─────────────────

println("\n── Markov model ───────────────────────────────────────────────────────")
markov_raw  = joinpath(RAW,  "markov")
markov_proc = joinpath(PROC, "markov")

for N in (1_000, 10_000), d in (0.0, 0.5, 0.9)
    fname = "neutral_markov_N$(N)_d$(d).jls"
    raw_path = joinpath(markov_raw, fname)
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, markov_proc, "markov N=$N d=$d")
end

# ── Deterministic files (2): N ∈ {1000, 10000}, d = 0.0 ──────────────────────

println("\n── Deterministic model ────────────────────────────────────────────────")
det_raw  = joinpath(RAW,  "deterministic")
det_proc = joinpath(PROC, "deterministic")

for N in (1_024, 16_384)
    fname = "neutral_deterministic_N$(N).jls"
    raw_path = joinpath(det_raw, fname)
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, det_proc, "deterministic N=$N")
end

println("\nAll done.")
