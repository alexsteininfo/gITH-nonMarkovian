using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Serialization
using Statistics
using Printf

const HELPERS = joinpath(dirname(dirname(@__DIR__)), "helpers")
include(joinpath(HELPERS, "types_selection1.jl"))

const RAW = joinpath(@__DIR__, "..", "..", "..", "data", "raw", "selection_1")

# (subdirectory, expected file count)
const EXPECTED = [("gamma", 120), ("markov", 120), ("deterministic", 40)]

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

        for x in sims
            isnothing(x.tree_root) &&
                push!(problems, "$sub/$f: nothing tree_root")
            x.injection.N_at_inject == x.params.N_critic + 1 ||
                push!(problems, "$sub/$f: N_at_inject $(x.injection.N_at_inject) ≠ N_critic+1 $(x.params.N_critic + 1)")
            x.injection.driver_clone_size > 0 ||
                push!(problems, "$sub/$f: driver clone size 0")
            x.params.nu == 2.0 ||
                push!(problems, "$sub/$f: nu $(x.params.nu) ≠ 2.0")
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

println()
println("total accepted simulations: $total_sims (expected 14000)")
if isempty(problems)
    println("INVENTORY OK")
else
    println("INVENTORY PROBLEMS ($(length(problems))):")
    foreach(p -> println("  - $p"), unique(problems))
    exit(1)
end
