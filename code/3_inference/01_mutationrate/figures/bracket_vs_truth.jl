# Figure 1: mutation_rate_bracket lower/upper endpoints vs. truth (m = 2).
#
# Reads every data/inference/mutationrate/neutral/full_trees/*.jls shard,
# groups by timing model (deterministic / gamma / markov), and draws a
# three-panel boxplot of lower/upper bracket endpoints with a horizontal
# reference line at M_TRUE = 2.0.
#
# Inputs:  data/inference/mutationrate/neutral/full_trees/*.jls
# Outputs: figures/3_inference/mutationrate/bracket_vs_truth.png
#
# Usage:
#   julia --project=. code/3_inference/01_mutationrate/figures/bracket_vs_truth.jl

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
using Pkg; Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics  # required to resolve SimParams / GrowthSimResult types in .jls files
using EvoTracer             # required to resolve StatTable / StatValue types in .jls files
using Serialization, Statistics
using CairoMakie
include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "plotting_functions.jl"))

const IN_DIR  = joinpath(DATA, "inference", "mutationrate", "neutral", "full_trees")
const OUT_DIR = joinpath(FIGURES, "3_inference", "mutationrate")
mkpath(OUT_DIR)
const M_TRUE = 2.0

# Parse filename stem → (timing, N, d, k). Deterministic files have no _k field.
function parse_stem(stem::String)
    s = replace(stem, r"^neutral_" => "")
    parts = split(s, "_")
    timing = parts[1]
    N = parse(Int,     replace(parts[2], "N" => ""))
    # Deterministic filenames omit _d and _k (d=0, k=Inf implicit).
    d = length(parts) >= 3 && startswith(parts[3], "d") ?
        parse(Float64, replace(parts[3], "d" => "")) : 0.0
    k = length(parts) >= 4 && startswith(parts[4], "k") ?
        parse(Float64, replace(parts[4], "k" => "")) : Inf
    return (; timing, N, d, k)
end

function collect_data()
    rows = Any[]
    isdir(IN_DIR) || return rows
    for f in sort(readdir(IN_DIR))
        endswith(f, ".jls") || continue
        stem = replace(f, ".jls" => "")
        pf = try parse_stem(stem) catch; continue end
        data = try deserialize(joinpath(IN_DIR, f)) catch e; @info "skip $f: $e"; continue end
        lowers, uppers = Float64[], Float64[]
        for d in data
            isnothing(d) && continue
            try
                tbl = d[:mutations][:mutation_rate_bracket]
                tbl isa Exception && continue
                bd = tbl.data.bound
                vl = tbl.data.value
                li = findfirst(==("lower"), bd)
                ui = findfirst(==("upper"), bd)
                isnothing(li) || isnothing(ui) && continue
                v_lo = vl[li]; v_up = vl[ui]
                isfinite(v_lo) && push!(lowers, v_lo)
                isfinite(v_up) && push!(uppers, v_up)
            catch e
                @info "skip sim in $stem: $e"
            end
        end
        push!(rows, (pf..., lowers=lowers, uppers=uppers))
    end
    return rows
end

function main()
    rows = collect_data()
    if isempty(rows)
        @info "No data found in $IN_DIR — producing empty figure."
        fig = Figure(size=(900, 350))
        Label(fig[1,1], "No data (IN_DIR empty)"; fontsize=FS_LABEL)
        save(joinpath(OUT_DIR, "bracket_vs_truth.png"), fig)
        println("wrote (empty) ", joinpath(OUT_DIR, "bracket_vs_truth.png"))
        return
    end

    timings = unique(r.timing for r in rows)
    sort!(timings)  # deterministic < gamma < markov alphabetically
    fig = Figure(size=(300 * length(timings), 400))
    Label(fig[0, :], "mutation_rate_bracket vs. truth  (m = $M_TRUE)";
          fontsize=FS_HEAD, halign=:center)

    timing_colors = Dict("deterministic" => COL_DET, "gamma" => COL_GAMMA, "markov" => COL_MARKOV)

    for (i, timing) in enumerate(timings)
        subs = filter(r -> r.timing == timing, rows)
        ax = Axis(fig[1, i];
                  title=timing, titlesize=FS_TITLE,
                  xlabel="shard", xlabelsize=FS_LABEL,
                  ylabel="m̂", ylabelsize=FS_LABEL,
                  xticklabelsize=FS_TICK, yticklabelsize=FS_TICK)
        hlines!(ax, [M_TRUE]; color=:red, linestyle=:dash, linewidth=2, label="m = $M_TRUE")

        for (j, r) in enumerate(subs)
            if !isempty(r.lowers)
                boxplot!(ax, fill(j - 0.2, length(r.lowers)), r.lowers;
                         width=0.35, color=(:steelblue, 0.65), strokewidth=1)
            end
            if !isempty(r.uppers)
                boxplot!(ax, fill(j + 0.2, length(r.uppers)), r.uppers;
                         width=0.35, color=(:darkorange, 0.65), strokewidth=1)
            end
        end

        xlabels = ["N=$(r.N)\nd=$(r.d)" for r in subs]
        ax.xticks = (collect(1:length(subs)), xlabels)
    end

    # Manual legend entries
    elem_lo  = PolyElement(color=(:steelblue,   0.65))
    elem_up  = PolyElement(color=(:darkorange,  0.65))
    elem_ref = LineElement(color=:red, linestyle=:dash)
    Legend(fig[2, :], [elem_lo, elem_up, elem_ref],
           ["lower endpoint", "upper endpoint", "truth m=$(M_TRUE)"];
           orientation=:horizontal, framevisible=false, labelsize=FS_LEGEND)

    out = joinpath(OUT_DIR, "bracket_vs_truth.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
