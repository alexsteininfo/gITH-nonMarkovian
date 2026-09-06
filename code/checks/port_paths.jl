# port_paths.jl — proves that the data and plot paths in a ported script resolve
# to exactly the same absolute locations as in its analysis/ original.
#
#   usage: julia --project=. code/checks/port_paths.jl <new-file> <old-file>
#
# This is the check that protects 23 GB. `port_equivalence.sh` deliberately
# strips the `const RAW|PROC|OUT|OUTDIR` lines, so it is blind to a path that
# changed meaning; this script evaluates those lines on both sides and compares.
#
# Method: find each `const NAME = joinpath(` declaration, extract the FULL
# joinpath(...) call by counting parens (some real constants — RAW_SHARD /
# SUB_SHARD in check_package_equivalence.jl — wrap their arguments onto a
# second physical line, so a single-line-only extraction mis-parses or
# crashes on exactly the constants Step 0 exists to cover), textually bind the
# names it depends on (@__DIR__ on the old side; DATA/PLOTS/ROOT on the new
# side), evaluate, and compare normalised absolute paths.
#
# TEMPORARY: delete this together with analysis/.

const REPO = dirname(dirname(@__DIR__))

const PATH_CONST = r"^[ \t]*const[ \t]+([A-Z_][A-Z_0-9]*)[ \t]*=[ \t]*(?=joinpath\()"m
const SKIP_NAMES = Set(["ROOT", "DATA", "PLOTS", "HELPERS"])

"""Starting at `start` (the index of the "j" in "joinpath("), return the
substring through the matching closing paren, however many lines it spans."""
function extract_joinpath(text::String, start::Int)
    depth = 0
    i = start
    while i <= ncodeunits(text)
        c = text[i]
        if c == '('
            depth += 1
        elseif c == ')'
            depth -= 1
            depth == 0 && return text[start:i]
        end
        i = nextind(text, i)
    end
    error("unbalanced joinpath( starting at byte $start")
end

# Some joinpath consts (RAW_SHARD / SUB_SHARD in check_package_equivalence.jl)
# interpolate another same-file constant, e.g. `joinpath(ROOT, ..., "$(STEM).jls")`.
# Bind every simple (non-joinpath) `const NAME = <literal>` elsewhere in the file
# as a plain variable first, so that name resolves when the joinpath expression
# is eval'd. Best-effort per line — a const this can't evaluate standalone is
# simply left unbound, and only matters if some joinpath expression actually
# needs it (in which case that expression's own eval fails loudly).
const SIMPLE_CONST = r"^[ \t]*const[ \t]+([A-Z_][A-Z_0-9]*)[ \t]*=[ \t]*(?!joinpath\()([^\n]+)$"m

function bind_simple_consts!(text::String)
    for m in eachmatch(SIMPLE_CONST, text)
        name, rhs = m.captures[1], m.captures[2]
        name in SKIP_NAMES && continue
        try
            Core.eval(Main, Meta.parse("$name = $rhs"))
        catch
        end
    end
end

"""Return Dict(name => resolved absolute path) for one file."""
function resolve_paths(file::String)
    dir  = dirname(abspath(file))
    text = read(file, String)
    bind_simple_consts!(text)
    out  = Dict{String,String}()
    for m in eachmatch(PATH_CONST, text)
        name  = m.captures[1]
        name in SKIP_NAMES && continue
        start = m.offset + ncodeunits(m.match)   # the lookahead consumed nothing
        expr  = extract_joinpath(text, start)
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
