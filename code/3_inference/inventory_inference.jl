# Sanity: does data/inference/ contain every file the parameter grid predicts?
# Runs in seconds; catches shard-naming drift and silently skipped writes.

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
using Pkg; Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

const STAGES     = ("mutationrate", "selection")
const SCENARIOS  = ("neutral", "selection_1", "selection_2")
const TREE_KINDS = ("full_trees", "subsampled_trees")

# Expected files are exactly the source files renamed under data/inference/…
function expected_for(scenario::String, tree_kind::String)
    src_root = tree_kind == "full_trees" ?
        joinpath(DATA, "raw",             scenario) :
        joinpath(DATA, "raw_subsampled",  scenario)
    files = String[]
    for (root, _, fs) in walkdir(src_root)
        for f in fs
            endswith(f, ".jls") || continue
            push!(files, f)
        end
    end
    return files
end

function main()
    total_missing = 0
    for stage in STAGES, scenario in SCENARIOS, tree_kind in TREE_KINDS
        expected  = Set(expected_for(scenario, tree_kind))
        out_dir   = joinpath(DATA, "inference", stage, scenario, tree_kind)
        present   = isdir(out_dir) ? Set(readdir(out_dir)) : Set{String}()
        missing_f = setdiff(expected, present)
        printstyled("  $stage/$scenario/$tree_kind: "; bold = true)
        println("$(length(present)) present / $(length(expected)) expected " *
                "($(length(missing_f)) missing)")
        total_missing += length(missing_f)
        if !isempty(missing_f) && length(missing_f) <= 10
            for m in sort(collect(missing_f))
                println("    missing: $m")
            end
        end
    end
    total_missing > 0 && error("$(total_missing) expected output files missing.")
    println("Inventory OK.")
end

isinteractive() || main()
