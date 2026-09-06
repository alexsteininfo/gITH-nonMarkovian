# port_paths.jl — proves that the data and plot paths in a ported script resolve
# to exactly the same absolute locations as in its analysis/ original.
#
#   usage: julia --project=. code/checks/port_paths.jl <new-file> <old-file>
#
# This is the check that protects 23 GB. `port_equivalence.sh` deliberately
# strips the `const RAW|PROC|OUT|OUTDIR` lines, so it is blind to a path that
# changed meaning; this script evaluates those lines on both sides and compares.
#
# Method: extract each `const NAME = joinpath(...)` line, textually bind the
# names it depends on (@__DIR__ on the old side; DATA/PLOTS/ROOT on the new
# side), evaluate, and compare normalised absolute paths.
#
# TEMPORARY: delete this together with analysis/.

const REPO = dirname(dirname(@__DIR__))

const PATH_CONST = r"^\s*const\s+([A-Z_][A-Z_0-9]*)\s*=\s*(joinpath\(.*)$"
const SKIP_NAMES = Set(["ROOT", "DATA", "PLOTS", "HELPERS"])

"""Return Dict(name => resolved absolute path) for one file."""
function resolve_paths(file::String)
    dir  = dirname(abspath(file))
    out  = Dict{String,String}()
    for line in eachline(file)
        m = match(PATH_CONST, line)
        isnothing(m) && continue
        name, expr = m.captures[1], m.captures[2]
        name in SKIP_NAMES && continue
        # Bind the free names textually (word-boundary aware, so a name like
        # RAW_ROOT is not corrupted by the ROOT substitution), then evaluate
        # the joinpath call.
        expr = replace(expr,
            "@__DIR__"    => repr(dir),
            r"\bDATA\b"   => repr(joinpath(REPO, "data")),
            r"\bPLOTS\b"  => repr(joinpath(REPO, "plots")),
            r"\bROOT\b"   => repr(REPO))
        out[name] = normpath(abspath(eval(Meta.parse(expr))))
    end
    return out
end

function main()
    length(ARGS) == 2 || error("usage: port_paths.jl <new-file> <old-file>")
    new_paths = resolve_paths(ARGS[1])
    old_paths = resolve_paths(ARGS[2])

    if keys(new_paths) != keys(old_paths)
        println("PATH CONSTANTS DIFFER: $(ARGS[1])")
        println("  analysis/: ", sort(collect(keys(old_paths))))
        println("  code/    : ", sort(collect(keys(new_paths))))
        exit(1)
    end

    bad = false
    for name in sort(collect(keys(old_paths)))
        if new_paths[name] != old_paths[name]
            println("PATH MISMATCH in $(ARGS[1]) — $name")
            println("  analysis/: ", old_paths[name])
            println("  code/    : ", new_paths[name])
            bad = true
        end
    end
    bad && exit(1)

    if isempty(old_paths)
        println("no path constants: $(ARGS[1])")
    else
        println("paths match ($(length(old_paths))): $(ARGS[1])")
    end
end

main()
