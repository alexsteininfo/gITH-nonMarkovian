using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using AbstractTrees
using Serialization
using Random

include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "subsampling.jl"))

const RAW       = joinpath(DATA, "raw_subsampled",       "selection_1")
const PROC      = joinpath(DATA, "processed_subsampled", "selection_1")
const FULL_PROC = joinpath(DATA, "processed",            "selection_1")

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
        error("$label: missing $full_injection_path — run code/data_generation/processing/process_sel1.jl first")

    # The full-tree params array lives alongside injection, with "params" substituted
    # for the quantity name.
    full_params_path = replace(full_injection_path, "injection" => "params")
    isfile(full_params_path) ||
        error("$label: missing $full_params_path — run code/data_generation/processing/process_sel1.jl first")

    println("  Processing: $(basename(raw_path))")

    subs      = deserialize(raw_path)::Vector{SubsampleResult{Sel1Params}}
    injection = deserialize(full_injection_path)::Vector{Sel1Injection}

    length(injection) == length(subs) ||
        error("$label: $(length(subs)) subsamples vs $(length(injection)) full-tree " *
              "injection entries — the two stages are not co-indexed")

    # The length check above passes under any permutation of the same multiset.
    # Params equality is fully discriminating: `Sel1Params` carries the
    # per-simulation `(N_critic, rep)`, so two co-indexed arrays can only match
    # element-for-element if every position lines up. The realistic trigger is
    # `process_sel1.jl` being resumable itself: a re-simulated raw shard combined
    # with a selectively-deleted full-tree output on one side and not the other
    # would still pass the length check but fail this one.
    full_params = deserialize(full_params_path)::Vector{Sel1Params}
    [s.params for s in subs] == full_params ||
        error("$label: subsampled params do not match $(basename(full_params_path)) " *
              "element-for-element for shard $(basename(raw_path)) — the subsampling " *
              "and full-tree processing stages have come apart")

    all_params       = Sel1Params[]
    all_mut_per_cell = Vector{Int}[]
    all_sfs          = Vector{Int}[]
    all_leaf_depths  = Vector{Int}[]
    all_leaf_fitness = Vector{Float64}[]

    for s in subs
        root   = s.tree_root
        depths = leaf_depths(root)

        length(depths) == s.n ||
            error("$label: sim $(s.sim_index) has $(length(depths)) leaves, expected n = $(s.n)")

        push!(all_params,       s.params)
        push!(all_mut_per_cell, mutations_per_cell(root))
        push!(all_sfs,          sitefrequencyspectrum(root, s.n))
        push!(all_leaf_depths,  depths)
        push!(all_leaf_fitness, leaf_fitness(root))
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
        # Atomic: the `isfile`-based resumability check above must never see a
        # truncated file left by a kill mid-write.
        serialize_atomic(joinpath(outdir, stem * ".jls"), data)
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
