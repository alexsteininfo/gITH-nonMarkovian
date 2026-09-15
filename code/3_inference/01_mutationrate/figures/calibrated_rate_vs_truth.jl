# Figure 3: calibrated_tree_length_rate median per shard vs. truth (m = 2).
#
# Reads every data/inference/mutationrate/neutral/full_trees/*.jls shard,
# computes the per-shard median of calibrated_tree_length_rate, and draws
# a scatter plot against M_TRUE = 2.0 (i.e., a jitter of points grouped by
# gamma_shape k, shaped by d), with a y = x diagonal reference line.
#
# Inputs:  data/inference/mutationrate/neutral/full_trees/*.jls
# Outputs: figures/3_inference/mutationrate/calibrated_rate_vs_truth.png
#
# Usage:
#   julia --project=. code/3_inference/01_mutationrate/figures/calibrated_rate_vs_truth.jl

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

function parse_stem(stem::String)
    s = replace(stem, r"^neutral_" => "")
    parts = split(s, "_")
    timing = parts[1]
    N = parse(Int,     replace(parts[2], "N" => ""))
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
        vals = Float64[]
        for d in data
            isnothing(d) && continue
            try
                sv = d[:mutations][:calibrated_tree_length_rate]
                sv isa Exception && continue
                isnothing(sv) && continue
                v = sv.value
                isfinite(v) && push!(vals, v)
            catch e
                @info "skip sim in $stem: $e"
            end
        end
        isempty(vals) && continue
        push!(rows, (pf..., med=median(vals), iqr_lo=quantile(vals, 0.25), iqr_hi=quantile(vals, 0.75)))
    end
    return rows
end

function main()
    rows = collect_data()
    fig = Figure(size=(700, 600))
    Label(fig[0, 1:2], "calibrated_tree_length_rate median vs. truth  (m = $M_TRUE)";
          fontsize=FS_HEAD, halign=:center)
    ax = Axis(fig[1, 1];
              xlabel="truth (m = $M_TRUE)", xlabelsize=FS_LABEL,
              ylabel="m̂_cal (median ± IQR)", ylabelsize=FS_LABEL,
              xticklabelsize=FS_TICK, yticklabelsize=FS_TICK)

    if isempty(rows)
        @info "No data found in $IN_DIR — producing empty figure."
        text!(ax, 0.5, 0.5; text="No data", align=(:center, :center),
              space=:relative, fontsize=FS_LABEL)
        save(joinpath(OUT_DIR, "calibrated_rate_vs_truth.png"), fig)
        println("wrote (empty) ", joinpath(OUT_DIR, "calibrated_rate_vs_truth.png"))
        return
    end

    # x-axis range: reference line
    all_meds = [r.med for r in rows if isfinite(r.med)]
    all_iqr  = vcat([r.iqr_hi for r in rows if isfinite(r.iqr_hi)],
                    [r.iqr_lo for r in rows if isfinite(r.iqr_lo)])
    ymax = maximum(vcat(all_meds, all_iqr, [M_TRUE])) * 1.1
    ymin = max(0, minimum(vcat(all_meds, all_iqr, [0.0])) * 0.9)
    xrange = LinRange(ymin, ymax, 50)
    lines!(ax, xrange, xrange; color=:black, linestyle=:dash, linewidth=2, label="y = x")

    # Color by k; marker shape by d
    all_k = sort(unique(r.k for r in rows if isfinite(r.k)))
    all_d = sort(unique(r.d for r in rows))
    k_colors = Dict(k => Makie.wong_colors()[i % 7 + 1] for (i, k) in enumerate(all_k))
    # push NaN k as a grey
    k_colors[NaN] = RGBAf(0.5, 0.5, 0.5, 1.0)
    d_markers = Dict(d => [:circle, :rect, :diamond, :utriangle, :dtriangle][i % 5 + 1]
                     for (i, d) in enumerate(all_d))

    for r in rows
        isfinite(r.med) || continue
        c  = get(k_colors, r.k, RGBAf(0.5, 0.5, 0.5, 1.0))
        mk = get(d_markers, r.d, :circle)
        errorbars!(ax, [M_TRUE], [r.med], [r.med - r.iqr_lo], [r.iqr_hi - r.med];
                   color=(c, 0.5), whiskerwidth=8)
        scatter!(ax, [M_TRUE], [r.med]; color=c, marker=mk, markersize=12, strokewidth=1)
    end

    # Legend: k (colors) + d (shapes)
    k_elems  = [MarkerElement(color=k_colors[k], marker=:circle, markersize=12) for k in all_k]
    k_labels = ["k=$(isnan(k) ? "det" : Int(k))" for k in all_k]
    d_elems  = [MarkerElement(color=RGBAf(0.5, 0.5, 0.5, 1.0), marker=d_markers[d], markersize=12) for d in all_d]
    d_labels = ["d=$(d)" for d in all_d]
    ref_elem = [LineElement(color=:black, linestyle=:dash)]

    Legend(fig[1, 2],
           vcat(k_elems, d_elems, ref_elem),
           vcat(k_labels, d_labels, ["y = x"]);
           framevisible=false, labelsize=FS_LEGEND)

    out = joinpath(OUT_DIR, "calibrated_rate_vs_truth.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
