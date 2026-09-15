# Figure 8: sibling_contrast_null.png
#
# Density plot of sibling contrast values pooled across neutral selection shards,
# one panel per timing model. The distribution should be symmetric around zero
# under the null. Heavier tails expected for Markov than Deterministic.
#
# Inputs:  data/inference/selection/neutral/full_trees/*.jls
# Outputs: figures/3_inference/selection/sibling_contrast_null.png
#
# Usage:
#   julia --project=. code/3_inference/02_selection/figures/sibling_contrast_null.jl

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
using Pkg; Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using EvoTracer
using Serialization, Statistics
using CairoMakie
include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "plotting_functions.jl"))

const IN_DIR  = joinpath(DATA, "inference", "selection", "neutral", "full_trees")
const OUT_DIR = joinpath(FIGURES, "3_inference", "selection")
mkpath(OUT_DIR)

function parse_stem_neutral(stem::String)
    s = replace(stem, r"^neutral_" => "")
    parts = split(s, "_")
    timing = parts[1]
    N = parse(Int, replace(parts[2], "N" => ""))
    d = length(parts) >= 3 && startswith(parts[3], "d") ?
        parse(Float64, replace(parts[3], "d" => "")) : 0.0
    k = length(parts) >= 4 && startswith(parts[4], "k") ?
        parse(Float64, replace(parts[4], "k" => "")) : Inf
    return (; timing, N, d, k)
end

function collect_data()
    by_timing = Dict{String, Vector{Float64}}()
    isdir(IN_DIR) || return by_timing
    for f in sort(readdir(IN_DIR))
        endswith(f, ".jls") || continue
        stem = replace(f, ".jls" => "")
        pf = try parse_stem_neutral(stem) catch; continue end
        timing = pf.timing
        data = try deserialize(joinpath(IN_DIR, f)) catch e; @info "skip $f: $e"; continue end
        vals = get!(by_timing, timing, Float64[])
        for d in data
            isnothing(d) && continue
            try
                tbl = d[:sibling_contrasts]
                (tbl isa Exception || isnothing(tbl)) && continue
                contrasts = tbl[:contrast]
                for v in contrasts
                    fv = Float64(v)
                    isfinite(fv) && push!(vals, fv)
                end
            catch e
                @info "skip sim in $stem: $e"
            end
        end
    end
    return by_timing
end

function main()
    by_timing = collect_data()
    if isempty(by_timing)
        @info "No data found in $IN_DIR — producing empty figure."
        fig = Figure(size=(700, 480))
        Label(fig[1, 1], "No data (IN_DIR empty)"; fontsize=FS_LABEL)
        out = joinpath(OUT_DIR, "sibling_contrast_null.png")
        save(out, fig)
        println("wrote (empty) ", out)
        return
    end

    timing_order  = ["deterministic", "gamma", "markov"]
    timing_colors = Dict("deterministic" => COL_DET, "gamma" => COL_GAMMA, "markov" => COL_MARKOV)
    timing_labels = Dict("deterministic" => "Deterministic", "gamma" => "Gamma (k=5)", "markov" => "Markov")

    fig = Figure(size=(700, 480))
    Label(fig[0, 1], "Sibling contrast null distribution (neutral trees)";
          fontsize=FS_HEAD, halign=:center)
    ax = Axis(fig[1, 1];
              xlabel="contrast δ_tot", xlabelsize=FS_LABEL,
              ylabel="density",        ylabelsize=FS_LABEL,
              xticklabelsize=FS_TICK, yticklabelsize=FS_TICK)

    legend_elements = []
    legend_labels   = String[]

    for timing in timing_order
        vals = get(by_timing, timing, Float64[])
        length(vals) < 2 && continue
        c = timing_colors[timing]
        density!(ax, vals; color=(c, 0.35), strokecolor=c, strokewidth=2)
        push!(legend_elements, PolyElement(color=(c, 0.5), strokecolor=c, strokewidth=1))
        push!(legend_labels, timing_labels[timing])
    end

    vlines!(ax, [0.0]; color=:black, linestyle=:dash, linewidth=2)
    push!(legend_elements, LineElement(color=:black, linestyle=:dash, linewidth=2))
    push!(legend_labels, "null (δ = 0)")

    if !isempty(legend_elements)
        Legend(fig[1, 2], legend_elements, legend_labels;
               framevisible=false, labelsize=FS_LEGEND)
    end

    out = joinpath(OUT_DIR, "sibling_contrast_null.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
