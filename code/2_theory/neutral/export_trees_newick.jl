using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using AbstractTrees
using Serialization
using Statistics

include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "newick.jl"))

# ── What gets exported ────────────────────────────────────────────────────────
#
# Four figures, each a 2x3 grid: columns are the timing models, rows are d = 0 / 0.9.
# The four form a 2x2 in (tips drawn) x (population sampled from) — see
# docs/superpowers/specs/2026-09-08-lineage-tree-figures-design.md.
#
# The deterministic sweep uses powers of two (1024, 16384) and has no d > 0 runs at
# all, so its d = 0.9 panel repeats its d = 0 panel and is flagged `repeated`.

const MODELS   = ["deterministic", "gamma", "markov"]
const D_VALUES = [0.0, 0.9]
const K        = 5.0

const FIGURE_KEYS = ["full_N1000", "n100_N1000", "n1000_N10000", "n100_N10000"]

# Hand-picked overrides, keyed (family, model, d). Empty by default: the median rule
# below chooses, so the figure is reproducible without a magic number.
const SIM_OVERRIDE = Dict{Tuple{Symbol,String,Float64},Int}()

const OUT = joinpath(DATA, "newick", "neutral")

# ── Source file names ─────────────────────────────────────────────────────────

function raw_path(model, N, d)
    model == "deterministic" && return joinpath(DATA, "raw", "neutral", "deterministic",
                                                "neutral_deterministic_N$(N).jls")
    model == "gamma" && return joinpath(DATA, "raw", "neutral", "gamma",
                                        "neutral_gamma_N$(N)_d$(d)_k$(K).jls")
    return joinpath(DATA, "raw", "neutral", "markov", "neutral_markov_N$(N)_d$(d).jls")
end

function sub_path(model, N, d, n)
    model == "deterministic" && return joinpath(DATA, "raw_subsampled", "neutral", "deterministic",
                                                "neutral_deterministic_N$(N)_n$(n).jls")
    model == "gamma" && return joinpath(DATA, "raw_subsampled", "neutral", "gamma",
                                        "neutral_gamma_N$(N)_d$(d)_k$(K)_n$(n).jls")
    return joinpath(DATA, "raw_subsampled", "neutral", "markov",
                    "neutral_markov_N$(N)_d$(d)_n$(n).jls")
end

# deterministic has no d > 0 runs; every d reads the d = 0 files
source_d(model, d)  = model == "deterministic" ? 0.0 : d
is_repeat(model, d) = model == "deterministic" && d != 0.0

# ── Simulation choice ─────────────────────────────────────────────────────────

"""
Index of the entry whose `t_end` is nearest the median, so a panel shows a typical
run rather than an outlier. Ties go to the lower index, which keeps the choice
reproducible.
"""
function nearest_median(t_ends::Vector{Float64})
    m = median(t_ends)
    return argmin(abs.(t_ends .- m))
end

# ── Panel rows ────────────────────────────────────────────────────────────────

struct PanelRow
    figure::String
    model::String
    d::Float64
    sim_index::Int
    n_tips::Int
    N_pop::Int
    sampling_fraction::Float64
    t_end::Float64
    root_edge::Float64
    median_leaf_depth::Int
    repeated::Bool
    newick_file::String
    highlight_file::String
end

# t_end and root_edge carry the same 12 significant digits as the branch lengths in
# the Newick files. R divides every node time by t_end, so rounding it — 6 decimals
# was the first attempt — pushes the tips just past 1.0 and makes "the tips reach the
# present" false at the precision the verifier checks. sampling_fraction is cosmetic.
csv_line(r::PanelRow) = join((
    r.figure, r.model, r.d, r.sim_index, r.n_tips, r.N_pop,
    round(r.sampling_fraction; digits = 6), _nwk_len(r.t_end),
    _nwk_len(r.root_edge), r.median_leaf_depth, r.repeated,
    r.newick_file, r.highlight_file), ',')

root_edge(root::BinaryNode{NonMarkovCell}) =
    (c = root.left === nothing ? root.right : root.left;
     c === nothing ? 0.0 : c.data.birthtime - root.data.birthtime)

"""
Write one tree and return its panel row. `newick_file` is a file *name*, not a path;
`repeated` panels are given a name that already exists and are not rewritten.
"""
function emit_panel(figure, model, d, root, tmax, sim_index, N_pop, force)
    nwk  = "$(figure)_$(model)_d$(source_d(model, d)).nwk"
    path = joinpath(OUT, nwk)
    if force || !isfile(path)
        write_newick(path, root; tmax = tmax)
    end
    depths = leaf_depths(root)                      # measured pre-collapse, on the tree
    n_tips = length(depths)
    hl = figure == "full_N1000" ? "highlight_$(model)_d$(source_d(model, d)).txt" : ""
    return PanelRow(figure, model, d, sim_index, n_tips, N_pop, n_tips / N_pop,
                    tmax, root_edge(root), round(Int, median(depths)),
                    is_repeat(model, d), nwk, hl)
end

# ── Main ──────────────────────────────────────────────────────────────────────

function main(force::Bool)
    mkpath(OUT)
    rows = PanelRow[]

    for model in MODELS, d in D_VALUES
        sd = source_d(model, d)

        # ---- the N ~ 1000 family: full tree + its n ~ 100 subsample ----
        N_small = model == "deterministic" ? 1_024 : 1_000
        n_small = model == "deterministic" ? 102   : 100

        raw   = deserialize(raw_path(model, N_small, sd))::Vector{GrowthSimResult}
        keep  = [(i, r) for (i, r) in enumerate(raw) if r.tree_root !== nothing]
        tends = [max(tree_tmax(r.tree_root), r.trajectory[end].t) for (_, r) in keep]

        subs     = deserialize(sub_path(model, N_small, sd, n_small))
        by_index = Dict(s.sim_index => s for s in subs)
        avail    = [j for (j, (i, _)) in enumerate(keep) if haskey(by_index, i)]
        isempty(avail) && error("no simulation is present in both $(raw_path(model, N_small, sd)) " *
                                "and its n = $(n_small) subsample")

        j = get(SIM_OVERRIDE, (:small, model, sd)) do
            avail[nearest_median(tends[avail])]
        end
        sim_index, result = keep[j]
        tmax = tends[j]
        println("  small  $(rpad(model, 14)) d=$(sd)  sim_index=$(sim_index)  " *
                "t_end=$(round(tmax, digits = 3))")

        push!(rows, emit_panel("full_N1000", model, d, result.tree_root, tmax,
                               sim_index, N_small, force))
        push!(rows, emit_panel("n100_N1000", model, d, by_index[sim_index].tree_root, tmax,
                               sim_index, N_small, force))

        # highlight: the tips of the full tree that the n ~ 100 sample retains
        hlfile = joinpath(OUT, "highlight_$(model)_d$(sd).txt")
        if force || !isfile(hlfile)
            write_atomic(hlfile, join(by_index[sim_index].sampled_ids, '\n') * "\n")
        end

        raw = nothing; keep = nothing; subs = nothing; by_index = nothing
        GC.gc()

        # ---- the N ~ 10000 family: two subsamples, no raw read ----
        N_big = model == "deterministic" ? 16_384 : 10_000
        n_hi  = model == "deterministic" ? 1_638  : 1_000
        n_lo  = model == "deterministic" ? 164    : 100

        big_hi  = deserialize(sub_path(model, N_big, sd, n_hi))
        hi_ends = [tree_tmax(s.tree_root) for s in big_hi]
        jb = get(SIM_OVERRIDE, (:big, model, sd)) do
            nearest_median(hi_ends)
        end
        big_index = big_hi[jb].sim_index
        tmax_big  = hi_ends[jb]
        println("  big    $(rpad(model, 14)) d=$(sd)  sim_index=$(big_index)  " *
                "t_end=$(round(tmax_big, digits = 3))")

        push!(rows, emit_panel("n1000_N10000", model, d, big_hi[jb].tree_root, tmax_big,
                               big_index, N_big, force))
        big_hi = nothing; GC.gc()

        big_lo = deserialize(sub_path(model, N_big, sd, n_lo))
        k = findfirst(s -> s.sim_index == big_index, big_lo)
        k === nothing && error("sim_index $(big_index) is missing from " *
                               sub_path(model, N_big, sd, n_lo))
        push!(rows, emit_panel("n100_N10000", model, d, big_lo[k].tree_root,
                               tree_tmax(big_lo[k].tree_root), big_index, N_big, force))
        big_lo = nothing; GC.gc()
    end

    # panels.csv is the index and is always rewritten, so it can never describe a
    # stale set of files. Row order follows the figure list, then model, then d.
    order = Dict(f => i for (i, f) in enumerate(FIGURE_KEYS))
    sort!(rows, by = r -> (order[r.figure], findfirst(==(r.model), MODELS), r.d))

    header = "figure,model,d,sim_index,n_tips,N_pop,sampling_fraction," *
             "t_end,root_edge,median_leaf_depth,repeated,newick_file,highlight_file"
    write_atomic(joinpath(OUT, "panels.csv"),
                 join(vcat(header, csv_line.(rows)), '\n') * "\n")

    println("\nWrote $(length(rows)) panel rows over " *
            "$(length(unique(r.newick_file for r in rows))) trees → $(OUT)")
end

# Already-complete runs exit before reading ~630 MB of trees.
const FORCE = "--force" in ARGS
if !FORCE && isfile(joinpath(OUT, "panels.csv"))
    done_rows = readlines(joinpath(OUT, "panels.csv"))[2:end]
    if !isempty(done_rows) && all(r -> isfile(joinpath(OUT, split(r, ',')[12])), done_rows)
        println("Already exported ($(length(done_rows)) panels) — pass --force to redo.")
        exit(0)
    end
end

main(FORCE)
