using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Serialization
using Statistics
using Printf

const HELPERS = joinpath(dirname(dirname(@__DIR__)), "helpers")
include(joinpath(HELPERS, "types_selection2.jl"))

const RAW = joinpath(@__DIR__, "..", "..", "..", "data", "raw", "selection_2")

# (subdirectory, expected file count)
#   gamma  / markov: 2 N_target × 3 d × 4 s = 24
#   deterministic:   2 N_target × 4 s       =  8   (d = 0 only)
const EXPECTED = [("gamma", 24), ("markov", 24), ("deterministic", 8)]

const N_SIMS       = 200    # simulations per file
const NU           = 1.0    # mutations per daughter per division
const M_CAP        = 10.0   # fitness cap
const EFFECT_SHAPE = 2.0    # a of the effect distribution Gamma(a, 1/a)

"""
    expected_from_filename(sub, f) -> NamedTuple

Parse the design point encoded in a shard's filename, for cross-checking against the
`Sel2Params` actually recorded inside the shard. Catches the one class of bug the
field-by-field checks cannot: a mismatch between what a file is named and what it holds.

The three shapes in use, each mirroring its paired neutral filename:

  - `sel2_gamma_N<N>_d<d>_k<k>_s<s>_M<M>.jls`     → `:gamma`,         gamma_shape from the file
  - `sel2_markov_N<N>_d<d>_s<s>_M<M>.jls`         → `:markov`,        gamma_shape 1.0
  - `sel2_deterministic_N<N>_s<s>_M<M>.jls`       → `:deterministic`, d 0.0, gamma_shape Inf
"""
function expected_from_filename(sub::AbstractString, f::AbstractString)
    if sub == "gamma"
        m = match(r"^sel2_gamma_N(\d+)_d([\d.]+)_k([\d.]+)_s([\d.]+)_M([\d.]+)\.jls$", f)
        isnothing(m) && error("gamma filename doesn't match expected pattern: $f")
        return (model = :gamma, N_target = parse(Int, m[1]), d = parse(Float64, m[2]),
                gamma_shape = parse(Float64, m[3]), s = parse(Float64, m[4]),
                M = parse(Float64, m[5]))
    elseif sub == "markov"
        m = match(r"^sel2_markov_N(\d+)_d([\d.]+)_s([\d.]+)_M([\d.]+)\.jls$", f)
        isnothing(m) && error("markov filename doesn't match expected pattern: $f")
        return (model = :markov, N_target = parse(Int, m[1]), d = parse(Float64, m[2]),
                gamma_shape = 1.0, s = parse(Float64, m[3]), M = parse(Float64, m[4]))
    elseif sub == "deterministic"
        m = match(r"^sel2_deterministic_N(\d+)_s([\d.]+)_M([\d.]+)\.jls$", f)
        isnothing(m) && error("deterministic filename doesn't match expected pattern: $f")
        return (model = :deterministic, N_target = parse(Int, m[1]), d = 0.0,
                gamma_shape = Inf, s = parse(Float64, m[2]), M = parse(Float64, m[3]))
    else
        error("unknown subdirectory: $sub")
    end
end

"""
    audit() -> (problems, total_sims, no_tree)

Body of the inventory sweep, wrapped in a function rather than run at top level: Julia's
soft-scope rules make a top-level `for` loop's `total_sims += ...` rebind a fresh local
on every iteration instead of accumulating, which crashes with `UndefVarError` before a
single row prints. A function body doesn't have that problem, and it's type-stable besides.

`nothing` tree roots are counted and reported rather than treated as problems: scenario 2
runs with `restart_on_extinction = true` and `run_sel2_once` documents that `getsingleroot`
may legitimately find no unique root, exactly as in the neutral runs, where
`process_neutral.jl` warns and skips such entries.
"""
function audit()
    problems   = String[]
    total_sims = 0
    no_tree    = 0

    for (sub, n_expected) in EXPECTED
        dir = joinpath(RAW, sub)
        if !isdir(dir)
            push!(problems, "missing directory $dir")
            continue
        end
        files = sort(filter(f -> endswith(f, ".jls"), readdir(dir)))
        length(files) == n_expected ||
            push!(problems, "$sub: $(length(files)) files, expected $n_expected")

        for f in files
            sims = deserialize(joinpath(dir, f))::Vector{Sel2SimResult}
            length(sims) == N_SIMS ||
                push!(problems, "$sub/$f: $(length(sims)) sims, expected $N_SIMS")

            expected  = expected_from_filename(sub, f)
            file_notree = 0

            for x in sims
                if isnothing(x.tree_root)
                    file_notree += 1
                    no_tree += 1
                end

                # Fixed background parameters.
                x.params.nu == NU ||
                    push!(problems, "$sub/$f: nu $(x.params.nu) ≠ $NU")
                x.params.effect_shape == EFFECT_SHAPE ||
                    push!(problems, "$sub/$f: effect_shape $(x.params.effect_shape) ≠ $EFFECT_SHAPE")
                x.params.b == 1.0 ||
                    push!(problems, "$sub/$f: b $(x.params.b) ≠ 1.0")

                # Cross-check the recorded params against what the filename promises.
                x.params.model == expected.model ||
                    push!(problems, "$sub/$f: model $(x.params.model) ≠ filename model $(expected.model)")
                x.params.N_target == expected.N_target ||
                    push!(problems, "$sub/$f: N_target $(x.params.N_target) ≠ filename N_target $(expected.N_target)")
                x.params.d == expected.d ||
                    push!(problems, "$sub/$f: d $(x.params.d) ≠ filename d $(expected.d)")
                x.params.s == expected.s ||
                    push!(problems, "$sub/$f: s $(x.params.s) ≠ filename s $(expected.s)")
                x.params.M == expected.M ||
                    push!(problems, "$sub/$f: M $(x.params.M) ≠ filename M $(expected.M)")
                x.params.gamma_shape == expected.gamma_shape ||
                    push!(problems, "$sub/$f: gamma_shape $(x.params.gamma_shape) ≠ filename gamma_shape $(expected.gamma_shape)")

                # Extinction/restart bookkeeping. With no death nothing can go extinct.
                x.n_restarts >= 0 ||
                    push!(problems, "$sub/$f: negative n_restarts $(x.n_restarts)")
                if x.params.d == 0.0 && x.n_restarts != 0
                    push!(problems, "$sub/$f: n_restarts $(x.n_restarts) ≠ 0 at d = 0")
                end

                # Fitness must stay in [1, M]: the update rule is non-decreasing from
                # f = 1 and the outer min caps it at M.
                if isempty(x.trajectory)
                    push!(problems, "$sub/$f: empty trajectory")
                else
                    fbar = last(x.trajectory).mean_fitness
                    (1.0 <= fbar <= x.params.M) ||
                        push!(problems, "$sub/$f: final mean fitness $fbar outside [1, $(x.params.M)]")
                end
            end

            # Every replicate present exactly once, numbered 1:N_SIMS.
            reps = sort([x.params.rep for x in sims])
            reps == collect(1:N_SIMS) ||
                push!(problems, "$sub/$f: reps are not exactly 1:$N_SIMS")
            allunique(x.params.seed for x in sims) ||
                push!(problems, "$sub/$f: duplicate seeds across replicates")

            total_sims += length(sims)

            rst  = mean(x.n_restarts for x in sims)
            fbar = mean(last(x.trajectory).mean_fitness
                        for x in sims if !isempty(x.trajectory))
            @printf("%-14s %-52s restarts %7.2f   mean fitness %6.2f   no-tree %3d\n",
                    sub, f, rst, fbar, file_notree)
        end
    end

    return problems, total_sims, no_tree
end

problems, total_sims, no_tree = audit()

const N_EXPECTED_TOTAL = sum(n for (_, n) in EXPECTED) * N_SIMS

println()
println("total simulations: $total_sims (expected $N_EXPECTED_TOTAL)")
println("simulations with no unique tree root: $no_tree (expected: some, skipped downstream)")
if isempty(problems)
    println("INVENTORY OK")
else
    unique_problems = unique(problems)
    println("INVENTORY PROBLEMS ($(length(unique_problems))):")
    foreach(p -> println("  - $p"), unique_problems)
    exit(1)
end
