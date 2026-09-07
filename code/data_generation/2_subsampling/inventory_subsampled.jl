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
include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "types_selection2.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "subsampling.jl"))

# Cross-check what is on disk against what the grids say should be there, in the
# style of `analysis/sim_runs/selection_1/inventory_sel1.jl`. Every expected
# filename is rebuilt from the same interpolation the writing scripts use, so a
# naming drift shows up as a missing file rather than as a silently ignored shard.

const S1 = collect(0.1:0.1:2.0)
const S2 = [0.05, 0.10, 0.15, 0.20]
const M  = 10.0

# Quick (default): the missing-file scan plus a 21-shard sample deep check, seconds.
# Full sweep (--deep): also deserializes every file under both data trees — several
# minutes over 4.75 GB. Printed up front so a reader of the output can never mistake
# one mode for the other.
const DEEP_MODE = "--deep" in ARGS
println(DEEP_MODE ?
    "Mode: FULL SWEEP (--deep) — every raw_subsampled and processed_subsampled file will be deserialized" :
    "Mode: quick (default) — missing-file scan + sampled deep check only; pass --deep for a full sweep")

# (scenario, model, full_stem, N_target)
function expected_stems()
    out = Tuple{String, String, String, Int}[]
    for N in (1_000, 10_000), d in (0.0, 0.5, 0.9)
        push!(out, ("neutral", "gamma",  "neutral_gamma_N$(N)_d$(d)_k5.0", N))
        push!(out, ("neutral", "markov", "neutral_markov_N$(N)_d$(d)",     N))
        for s in S1
            push!(out, ("selection_1", "gamma",  "sel1_gamma_N$(N)_d$(d)_k5.0_s$(s)", N))
            push!(out, ("selection_1", "markov", "sel1_markov_N$(N)_d$(d)_s$(s)",     N))
        end
        for s in S2
            push!(out, ("selection_2", "gamma",  "sel2_gamma_N$(N)_d$(d)_k5.0_s$(s)_M$(M)", N))
            push!(out, ("selection_2", "markov", "sel2_markov_N$(N)_d$(d)_s$(s)_M$(M)",     N))
        end
    end
    for N in (1_024, 16_384)
        push!(out, ("neutral", "deterministic", "neutral_deterministic_N$(N)", N))
        for s in S1
            push!(out, ("selection_1", "deterministic", "sel1_deterministic_N$(N)_s$(s)", N))
        end
        for s in S2
            push!(out, ("selection_2", "deterministic", "sel2_deterministic_N$(N)_s$(s)_M$(M)", N))
        end
    end
    return out
end

const QUANTITIES = Dict(
    "neutral"     => ("params", "mut_per_cell", "sfs", "leaf_depths"),
    "selection_1" => ("params", "mut_per_cell", "sfs", "leaf_depths", "leaf_fitness", "injection"),
    "selection_2" => ("params", "mut_per_cell", "sfs", "leaf_depths", "leaf_fitness", "n_restarts"),
)

missing_raw  = String[]
missing_proc = String[]
n_raw = 0
n_proc = 0

for (scenario, model, stem, N) in expected_stems()
    global n_raw, n_proc
    for n in sample_sizes(N)
        name = "$(stem)_n$(n).jls"

        raw = joinpath(ROOT, "data", "raw_subsampled", scenario, model, name)
        isfile(raw) ? (n_raw += 1) : push!(missing_raw, raw)

        for q in QUANTITIES[scenario]
            p = joinpath(ROOT, "data", "processed_subsampled", scenario, model, q, name)
            isfile(p) ? (n_proc += 1) : push!(missing_proc, p)
        end
    end
end

println("raw_subsampled       : $n_raw files present, $(length(missing_raw)) missing")
println("processed_subsampled : $n_proc files present, $(length(missing_proc)) missing")

for p in first(missing_raw,  20); println("  MISSING raw:  ", relpath(p, ROOT)); end
for p in first(missing_proc, 20); println("  MISSING proc: ", relpath(p, ROOT)); end

# Deep-check one raw shard per (scenario, model, n) combination — 21 of 525, first
# simulation only, and no processed file is opened here (that's what --deep below is
# for). For the shard it does open: the first entry's recorded n and N_full match
# the filename, and its tree's leaves are exactly the recorded sampled_ids.
println("\n── Deep check (quick sample) ─────────────────────────────────────────")
seen = Set{Tuple{String, String, Int}}()
problems = 0
for (scenario, model, stem, N) in expected_stems()
    global problems
    for n in sample_sizes(N)
        key = (scenario, model, n)
        key in seen && continue
        raw = joinpath(ROOT, "data", "raw_subsampled", scenario, model, "$(stem)_n$(n).jls")
        isfile(raw) || continue
        push!(seen, key)

        # A corrupt/truncated file must be tallied as a malformed shard, not abort
        # the whole run with a stack trace — this is the only place a genuinely
        # unreadable raw shard would otherwise go unnoticed in the default quick mode.
        try
            subs = deserialize(raw)
            r    = subs[1]
            ids  = sort([l.data.id for l in Leaves(r.tree_root)])
            ok   = r.n == n && r.N_full == N && ids == sort(r.sampled_ids) &&
                   length(r.sampled_ids) == n && length(unique(r.sampled_ids)) == n
            ok || (problems += 1)
            println(rpad("$scenario/$model n=$n", 40), ok ? " OK" : " PROBLEM",
                    "  sims=", length(subs), "  eltype=", eltype(subs))
        catch e
            problems += 1
            println(rpad("$scenario/$model n=$n", 40), " MALFORMED  ",
                    relpath(raw, ROOT), "  (", typeof(e), ")")
        end
    end
end

# ── Full sweep (opt-in via --deep) ──────────────────────────────────────────────
#
# The quick deep check above proves almost nothing about the bulk of the data: it
# opens 21 of 525 raw shards and 0 of 3108 processed files. This instead
# deserializes every file under both data trees and only checks that the read
# succeeds — proof against truncation/corruption, not a repeat of the quick check's
# structural validation, which already covers a sample.
function full_sweep()
    n_swept = 0
    n_unreadable = 0
    for dir in (joinpath(ROOT, "data", "raw_subsampled"),
                joinpath(ROOT, "data", "processed_subsampled"))
        for (dirpath, _, files) in walkdir(dir)
            for f in files
                endswith(f, ".jls") || continue
                path = joinpath(dirpath, f)
                n_swept += 1
                try
                    deserialize(path)
                catch e
                    n_unreadable += 1
                    println("  UNREADABLE: ", relpath(path, ROOT), "  (", typeof(e), ")")
                end
            end
        end
    end
    return n_swept, n_unreadable
end

n_unreadable = 0
if DEEP_MODE
    println("\n── Full sweep (--deep) ─────────────────────────────────────────────────")
    n_swept, n_unreadable = full_sweep()
    println("full sweep: $n_swept files read, $n_unreadable unreadable")
end

total_missing = length(missing_raw) + length(missing_proc)
problems += n_unreadable
println("\n$(total_missing) missing files, $problems malformed shards",
        DEEP_MODE ? " (includes $n_unreadable unreadable in the full sweep)" : "")
println(DEEP_MODE ? "Mode was: FULL SWEEP (--deep)" : "Mode was: quick (default)")
exit(total_missing == 0 && problems == 0 ? 0 : 1)
