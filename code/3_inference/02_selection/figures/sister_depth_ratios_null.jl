# Figure 9: sister_depth_ratios_null.png
#
# Density plot of sister depth ratios pooled across neutral MUTATION-RATE shards,
# one curve per timing model. Plotted on a log x-axis; distributions should be
# centred on 1 (log(ratio) = 0). Wider spread expected for Markov than Deterministic.
#
# NOTE: This reads mutation-rate outputs (data/inference/mutationrate/), not
# selection outputs, because sister_depth_ratios is computed as part of the
# mutation-rate suite.
#
# The relevant column from sister_depth_ratios is :rate_ratio (see ratecontrasts.jl).
#
# Inputs:  data/inference/mutationrate/neutral/full_trees/*.jls
# Outputs: figures/3_inference/selection/sister_depth_ratios_null.png
#
# Usage:
#   julia --project=. code/3_inference/02_selection/figures/sister_depth_ratios_null.jl

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

const IN_DIR  = joinpath(DATA, "inference", "mutationrate", "neutral", "full_trees")
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
                # mutation-rate output: d[:mutations][:sister_depth_ratios]
                # column is :rate_ratio (from ratecontrasts.jl docstring)
                mutations = d[:mutations]
                (mutations isa Exception || isnothing(mutations)) && continue
                tbl = mutations[:sister_depth_ratios]
                (tbl isa Exception || isnothing(tbl)) && continue
                # Try :rate_ratio first, fall back to :ratio for forward compatibility
                col_name = :rate_ratio in propertynames(tbl.data) ? :rate_ratio : :ratio
                ratios = tbl[col_name]
                for v in ratios
                    fv = Float64(v)
                    (isfinite(fv) && fv > 0) && push!(vals, fv)
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
        out = joinpath(OUT_DIR, "sister_depth_ratios_null.png")
        save(out, fig)
        println("wrote (empty) ", out)
        return
    end

    timing_order  = ["deterministic", "gamma", "markov"]
    timing_colors = Dict("deterministic" => COL_DET, "gamma" => COL_GAMMA, "markov" => COL_MARKOV)
    timing_labels = Dict("deterministic" => "Deterministic", "gamma" => "Gamma (k=5)", "markov" => "Markov")

    fig = Figure(size=(700, 480))
    Label(fig[0, 1], "Sister depth ratio null distribution (neutral trees, log₁₀ x)";
          fontsize=FS_HEAD, halign=:center)
    # Plot log10(rate_ratio) on a linear axis to avoid Makie's log-scale density issue.
    ax = Axis(fig[1, 1];
              xlabel="log₁₀(depth rate_ratio)", xlabelsize=FS_LABEL,
              ylabel="density",                   ylabelsize=FS_LABEL,
              xticklabelsize=FS_TICK, yticklabelsize=FS_TICK)

    legend_elements = []
    legend_labels   = String[]

    for timing in timing_order
        vals = get(by_timing, timing, Float64[])
        length(vals) < 2 && continue
        log_vals = log10.(filter(v -> v > 0, vals))
        length(log_vals) < 2 && continue
        c = timing_colors[timing]
        density!(ax, log_vals; color=(c, 0.35), strokecolor=c, strokewidth=2)
        push!(legend_elements, PolyElement(color=(c, 0.5), strokecolor=c, strokewidth=1))
        push!(legend_labels, timing_labels[timing])
    end

    vlines!(ax, [0.0]; color=:red, linestyle=:dash, linewidth=2)  # log10(1) = 0
    push!(legend_elements, LineElement(color=:red, linestyle=:dash, linewidth=2))
    push!(legend_labels, "null (ratio = 1, log = 0)")

    if !isempty(legend_elements)
        Legend(fig[1, 2], legend_elements, legend_labels;
               framevisible=false, labelsize=FS_LEGEND)
    end

    out = joinpath(OUT_DIR, "sister_depth_ratios_null.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
