using Pkg
Pkg.activate(dirname(dirname(@__DIR__)))

using MutationLoadDynamics
using AbstractTrees
using Serialization
using Random

const HELPERS = joinpath(dirname(@__DIR__), "helpers")

include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "tree_analysis.jl"))
include(joinpath(HELPERS, "subsampling.jl"))

const RAW       = joinpath(@__DIR__, "..", "..", "data", "raw_subsampled",       "selection_1")
const PROC      = joinpath(@__DIR__, "..", "..", "data", "processed_subsampled", "selection_1")
const FULL_PROC = joinpath(@__DIR__, "..", "..", "data", "processed",            "selection_1")

# Stage 2b for selection scenario 1, recomputed on a random subsample of n cells.
#
#   params       — Sel1Params of the source simulation
#   mut_per_cell — per-cell neutral burden of each sampled cell
#   sfs          — sfs[k] = mutations carried by exactly k of the n sampled cells
#   leaf_depths  — divisions from founder to each sampled cell
#   leaf_fitness — final fitness of each sampled cell, co-indexed with mut_per_cell
#   injection    — Sel1Injection of the source simulation
#
# `leaf_fitness` is stored here even though `process_sel1.jl` does not store it for
# full trees. Scenario 1's headline observable is `injection.driver_clone_size`, the
# number of cells with fitness > 1 at N_target; its sample analogue is the number of
# *sampled* cells with fitness > 1, and that cannot be recovered from any other
# stored quantity — the driver deliberately does not increment `mutations`, so the
# clone is invisible to `sfs` and `mut_per_cell`. The count itself is left to the
# plotting stage rather than baked in here.
#
# `injection` is copied from `data/processed/selection_1/`, not re-derived from
# `data/raw/`: those arrays were written by `process_sel1.jl` under the same
# nothing-tree filter in the same source order, so they are already
# element-for-element aligned, and re-deriving would mean a second full pass over
# 6.9 GB of trees. This stage therefore requires
# `analysis/processing/process_sel1.jl` to have run first.
#
# Each shard holds 50 sims: 10 log-spaced N_critic × 5 reps at one (model, N, d, s).

const SUBDIRS = ("params", "mut_per_cell", "sfs", "leaf_depths", "leaf_fitness", "injection")

function process_file(raw_path::String, proc_dir::String, full_injection_path::String,
                      label::String)
    stem = splitext(basename(raw_path))[1]

    if all(isfile(joinpath(proc_dir, sub, stem * ".jls")) for sub in SUBDIRS)
        println("  Skipping (already processed): $(basename(raw_path))")
        return
    end

    isfile(full_injection_path) ||
        error("$label: missing $full_injection_path — run analysis/processing/process_sel1.jl first")

    println("  Processing: $(basename(raw_path))")

    subs      = deserialize(raw_path)::Vector{SubsampleResult{Sel1Params}}
    injection = deserialize(full_injection_path)::Vector{Sel1Injection}

    length(injection) == length(subs) ||
        error("$label: $(length(subs)) subsamples vs $(length(injection)) full-tree " *
              "injection entries — the two stages are not co-indexed")

    all_params       = Sel1Params[]
    all_mut_per_cell = Vector{Int}[]
    all_sfs          = Vector{Int}[]
    all_leaf_depths  = Vector{Int}[]
    all_leaf_fitness = Vector{Float64}[]

    for s in subs
        root   = s.tree_root
        depths = compute_leaf_depths(root)

        length(depths) == s.n ||
            error("$label: sim $(s.sim_index) has $(length(depths)) leaves, expected n = $(s.n)")

        push!(all_params,       s.params)
        push!(all_mut_per_cell, compute_mut_per_cell(root))
        push!(all_sfs,          compute_sfs(root, s.n))
        push!(all_leaf_depths,  depths)
        push!(all_leaf_fitness, compute_leaf_fitness(root))
    end

    for (subdir, data) in (
        ("params",        all_params),
        ("mut_per_cell",  all_mut_per_cell),
        ("sfs",           all_sfs),
        ("leaf_depths",   all_leaf_depths),
        ("leaf_fitness",  all_leaf_fitness),
        ("injection",     injection),
    )
        outdir = joinpath(proc_dir, subdir)
        mkpath(outdir)
        serialize(joinpath(outdir, stem * ".jls"), data)
    end

    println("    → saved $(length(all_params)) results to $proc_dir")
end

const S_VALUES = collect(0.1:0.1:2.0)

println("\n── Gamma model ────────────────────────────────────────────────────────")
for N in (1_000, 10_000), d in (0.0, 0.5, 0.9), s in S_VALUES, n in sample_sizes(N)
    full_stem = "sel1_gamma_N$(N)_d$(d)_k5.0_s$(s)"
    raw_path  = joinpath(RAW, "gamma", "$(full_stem)_n$(n).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, joinpath(PROC, "gamma"),
                 joinpath(FULL_PROC, "gamma", "injection", full_stem * ".jls"),
                 "gamma N=$N d=$d s=$s n=$n")
end

println("\n── Markov model ───────────────────────────────────────────────────────")
for N in (1_000, 10_000), d in (0.0, 0.5, 0.9), s in S_VALUES, n in sample_sizes(N)
    full_stem = "sel1_markov_N$(N)_d$(d)_s$(s)"
    raw_path  = joinpath(RAW, "markov", "$(full_stem)_n$(n).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, joinpath(PROC, "markov"),
                 joinpath(FULL_PROC, "markov", "injection", full_stem * ".jls"),
                 "markov N=$N d=$d s=$s n=$n")
end

println("\n── Deterministic model ────────────────────────────────────────────────")
for N in (1_024, 16_384), s in S_VALUES, n in sample_sizes(N)
    full_stem = "sel1_deterministic_N$(N)_s$(s)"
    raw_path  = joinpath(RAW, "deterministic", "$(full_stem)_n$(n).jls")
    isfile(raw_path) || (@warn "Not found: $raw_path"; continue)
    process_file(raw_path, joinpath(PROC, "deterministic"),
                 joinpath(FULL_PROC, "deterministic", "injection", full_stem * ".jls"),
                 "deterministic N=$N s=$s n=$n")
end

println("\nAll done.")
