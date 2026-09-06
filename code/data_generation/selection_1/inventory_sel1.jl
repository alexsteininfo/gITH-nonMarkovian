using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using Serialization
using Statistics
using Printf

include(joinpath(HELPERS, "types_selection1.jl"))

const RAW = joinpath(DATA, "raw", "selection_1")

# (subdirectory, expected file count)
const EXPECTED = [("gamma", 120), ("markov", 120), ("deterministic", 40)]

"""
    expected_from_filename(sub, f) -> NamedTuple

Parse the design point encoded in a shard's filename, for cross-checking against the
`Sel1Params` actually recorded inside the shard. Catches exactly the filename / transposed
argument bugs the field-by-field checks below can't: a mismatch between what the file is
named and what it contains.

The three shapes in use:

  - `sel1_gamma_N<N>_d<d>_k<k>_s<s>.jls`         → model `:gamma`,         gamma_shape 5.0
  - `sel1_markov_N<N>_d<d>_s<s>.jls`              → model `:markov`,       gamma_shape 1.0
  - `sel1_deterministic_N<N>_s<s>.jls`            → model `:deterministic`, d 0.0, gamma_shape Inf
"""
function expected_from_filename(sub::AbstractString, f::AbstractString)
    if sub == "gamma"
        m = match(r"^sel1_gamma_N(\d+)_d([\d.]+)_k([\d.]+)_s([\d.]+)\.jls$", f)
        isnothing(m) && error("gamma filename doesn't match expected pattern: $f")
        return (model = :gamma, N_target = parse(Int, m[1]), d = parse(Float64, m[2]),
                gamma_shape = parse(Float64, m[3]), s = parse(Float64, m[4]))
    elseif sub == "markov"
        m = match(r"^sel1_markov_N(\d+)_d([\d.]+)_s([\d.]+)\.jls$", f)
        isnothing(m) && error("markov filename doesn't match expected pattern: $f")
        return (model = :markov, N_target = parse(Int, m[1]), d = parse(Float64, m[2]),
                gamma_shape = 1.0, s = parse(Float64, m[3]))
    elseif sub == "deterministic"
        m = match(r"^sel1_deterministic_N(\d+)_s([\d.]+)\.jls$", f)
        isnothing(m) && error("deterministic filename doesn't match expected pattern: $f")
        return (model = :deterministic, N_target = parse(Int, m[1]), d = 0.0,
                gamma_shape = Inf, s = parse(Float64, m[2]))
    else
        error("unknown subdirectory: $sub")
    end
end

"""
    audit() -> (problems, total_sims)

Body of the inventory sweep, wrapped in a function rather than run at top level: Julia's
soft-scope rules make a top-level `for` loop's `total_sims += ...` rebind a fresh local
on every iteration instead of accumulating, which used to crash the script with
`UndefVarError` before it printed a single row. A function body doesn't have that
problem, and it's type-stable besides.
"""
function audit()
    problems = String[]
    total_sims = 0

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
            sims = deserialize(joinpath(dir, f))::Vector{Sel1SimResult}
            length(sims) == 50 ||
                push!(problems, "$sub/$f: $(length(sims)) sims, expected 50")

            expected = expected_from_filename(sub, f)

            for x in sims
                isnothing(x.tree_root) &&
                    push!(problems, "$sub/$f: nothing tree_root")
                x.injection.N_at_inject == x.params.N_critic + 1 ||
                    push!(problems, "$sub/$f: N_at_inject $(x.injection.N_at_inject) ≠ N_critic+1 $(x.params.N_critic + 1)")
                x.injection.driver_clone_size > 0 ||
                    push!(problems, "$sub/$f: driver clone size 0")
                x.injection.driver_cell_id > 0 ||
                    push!(problems, "$sub/$f: driver_cell_id $(x.injection.driver_cell_id) not positive")
                x.params.nu == 2.0 ||
                    push!(problems, "$sub/$f: nu $(x.params.nu) ≠ 2.0")

                # Cross-check the recorded params against what the filename promises.
                x.params.model == expected.model ||
                    push!(problems, "$sub/$f: model $(x.params.model) ≠ filename model $(expected.model)")
                x.params.N_target == expected.N_target ||
                    push!(problems, "$sub/$f: N_target $(x.params.N_target) ≠ filename N_target $(expected.N_target)")
                x.params.d == expected.d ||
                    push!(problems, "$sub/$f: d $(x.params.d) ≠ filename d $(expected.d)")
                x.params.s == expected.s ||
                    push!(problems, "$sub/$f: s $(x.params.s) ≠ filename s $(expected.s)")
                x.params.gamma_shape == expected.gamma_shape ||
                    push!(problems, "$sub/$f: gamma_shape $(x.params.gamma_shape) ≠ filename gamma_shape $(expected.gamma_shape)")
            end

            # every (N_critic, rep) cell present exactly once
            cells = sort([(x.params.N_critic, x.params.rep) for x in sims])
            allunique(cells) || push!(problems, "$sub/$f: duplicate (N_critic, rep) cells")
            length(unique(x.params.N_critic for x in sims)) == 10 ||
                push!(problems, "$sub/$f: not 10 distinct N_critic values")

            total_sims += length(sims)

            att = mean(x.injection.n_attempts for x in sims)
            cln = mean(x.injection.driver_clone_size for x in sims)
            @printf("%-14s %-46s attempts %6.2f   mean clone %9.1f\n", sub, f, att, cln)
        end
    end

    return problems, total_sims
end

problems, total_sims = audit()

println()
println("total accepted simulations: $total_sims (expected 14000)")
if isempty(problems)
    println("INVENTORY OK")
else
    unique_problems = unique(problems)
    println("INVENTORY PROBLEMS ($(length(unique_problems))):")
    foreach(p -> println("  - $p"), unique_problems)
    exit(1)
end
