# Figure 2: cherry_contrast distribution, one density curve per timing model.
#
# Reads every data/inference/mutationrate/neutral/full_trees/*.jls shard,
# groups by timing model, and draws overlaid density curves of cherry_contrast
# with a vertical dashed line at x = 1 (the Poisson-clock null).
#
# Inputs:  data/inference/mutationrate/neutral/full_trees/*.jls
# Outputs: figures/3_inference/mutationrate/cherry_contrast_by_timing.png
#
# Usage:
#   julia --project=. code/3_inference/01_mutationrate/figures/cherry_contrast_by_timing.jl

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

function parse_timing(stem::String)
    m = match(r"^neutral_(deterministic|gamma|markov)_", stem)
    isnothing(m) ? "unknown" : m.captures[1]
end

function collect_data()
    by_timing = Dict{String, Vector{Float64}}()
    isdir(IN_DIR) || return by_timing
    for f in sort(readdir(IN_DIR))
        endswith(f, ".jls") || continue
        stem = replace(f, ".jls" => "")
        timing = parse_timing(stem)
        data = try deserialize(joinpath(IN_DIR, f)) catch e; @info "skip $f: $e"; continue end
        vals = get!(by_timing, timing, Float64[])
        for d in data
            isnothing(d) && continue
            try
                sv = d[:mutations][:cherry_contrast]
                sv isa Exception && continue
                v = sv.value
                isfinite(v) && push!(vals, v)
            catch e
                @info "skip sim in $stem: $e"
            end
        end
    end
    return by_timing
end

function main()
    by_timing = collect_data()
    fig = Figure(size=(700, 480))
    Label(fig[0, 1], "cherry_contrast by timing model";
          fontsize=FS_HEAD, halign=:center)
    ax = Axis(fig[1, 1];
              xlabel="cherry contrast ratio", xlabelsize=FS_LABEL,
              ylabel="density",               ylabelsize=FS_LABEL,
              xticklabelsize=FS_TICK, yticklabelsize=FS_TICK)

    timing_order  = ["deterministic", "gamma", "markov"]
    timing_colors = Dict("deterministic" => COL_DET, "gamma" => COL_GAMMA, "markov" => COL_MARKOV)
    timing_labels = Dict("deterministic" => "Deterministic", "gamma" => "Gamma (k=5)", "markov" => "Markov")

    legend_elements = []
    legend_labels   = String[]

    for timing in timing_order
        vals = get(by_timing, timing, Float64[])
        length(vals) < 2 && continue
        c = timing_colors[timing]
        density!(ax, vals; color=(c, 0.35), strokecolor=c, strokewidth=2)
        push!(legend_elements, PolyElement(color=(c, 0.5), strokecolor=c, strokewidth=1))
        push!(legend_labels,   timing_labels[timing])
    end

    vlines!(ax, [1.0]; color=:red, linestyle=:dash, linewidth=2, label="Poisson-clock null")
    push!(legend_elements, LineElement(color=:red, linestyle=:dash, linewidth=2))
    push!(legend_labels,   "null (ratio = 1)")

    if !isempty(legend_elements)
        Legend(fig[1, 2], legend_elements, legend_labels;
               framevisible=false, labelsize=FS_LEGEND)
    end

    if isempty(by_timing)
        @info "No data found in $IN_DIR — producing empty figure."
        text!(ax, 0.5, 0.5; text="No data", align=(:center, :center),
              space=:relative, fontsize=FS_LABEL)
    end

    out = joinpath(OUT_DIR, "cherry_contrast_by_timing.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
