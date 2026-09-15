# Figure 7: splitting_null_pvalue_distribution.png
#
# Histograms of splitting-null p-values pooled across neutral selection shards,
# one panel per timing model. A red dashed horizontal line marks density = 1,
# the expected uniform density under the null. Rows are split by :exact flag
# (binary exact null vs. polytomy provisional null) when both are present.
#
# Inputs:  data/inference/selection/neutral/full_trees/*.jls
# Outputs: figures/3_inference/selection/splitting_null_pvalue_distribution.png
#
# Usage:
#   julia --project=. code/3_inference/02_selection/figures/splitting_null_pvalue_distribution.jl

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
    # by_timing => (exact_pvals, simulated_pvals)
    by_timing = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
    isdir(IN_DIR) || return by_timing
    for f in sort(readdir(IN_DIR))
        endswith(f, ".jls") || continue
        stem = replace(f, ".jls" => "")
        pf = try parse_stem_neutral(stem) catch; continue end
        timing = pf.timing
        data = try deserialize(joinpath(IN_DIR, f)) catch e; @info "skip $f: $e"; continue end
        exact_v, sim_v = get!(by_timing, timing,
                               (Float64[], Float64[]))
        for d in data
            isnothing(d) && continue
            try
                tbl = d[:splitting_null_pvalues]
                (tbl isa Exception || isnothing(tbl)) && continue
                pvals  = tbl[:pvalue]
                exacts = tbl[:exact]
                for i in eachindex(pvals)
                    v = Float64(pvals[i])
                    isfinite(v) || continue
                    if Bool(exacts[i])
                        push!(exact_v, v)
                    else
                        push!(sim_v, v)
                    end
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
        fig = Figure(size=(900, 400))
        Label(fig[1, 1], "No data (IN_DIR empty)"; fontsize=FS_LABEL)
        out = joinpath(OUT_DIR, "splitting_null_pvalue_distribution.png")
        save(out, fig)
        println("wrote (empty) ", out)
        return
    end

    timing_order  = ["deterministic", "gamma", "markov"]
    timing_colors = Dict("deterministic" => COL_DET, "gamma" => COL_GAMMA, "markov" => COL_MARKOV)
    timing_labels = Dict("deterministic" => "Deterministic", "gamma" => "Gamma (k=5)", "markov" => "Markov")
    timings = filter(t -> haskey(by_timing, t), timing_order)

    # Determine if we have polytomy rows to show
    has_polytomy = any(t -> !isempty(by_timing[t][2]), timings)
    nrows = has_polytomy ? 2 : 1
    ncols = length(timings)

    fig = Figure(size=(350 * ncols, 350 * nrows + 60))
    Label(fig[0, :], "Splitting null p-value distribution (neutral trees)";
          fontsize=FS_HEAD, halign=:center)

    for (j, timing) in enumerate(timings)
        exact_v, sim_v = by_timing[timing]
        c = timing_colors[timing]
        lbl = timing_labels[timing]

        # Row 1: exact (binary) null
        ax1 = Axis(fig[1, j];
                   title="$lbl — binary (exact)",
                   titlesize=FS_TITLE,
                   xlabel="p-value", xlabelsize=FS_LABEL,
                   ylabel="density",  ylabelsize=FS_LABEL,
                   xticklabelsize=FS_TICK, yticklabelsize=FS_TICK,
                   limits=((0, 1), nothing))
        if !isempty(exact_v)
            hist!(ax1, exact_v; normalization=:pdf, bins=20,
                  color=(c, 0.5), strokewidth=1, strokecolor=c)
        else
            text!(ax1, 0.5, 0.5; text="no data", align=(:center, :center),
                  space=:relative, fontsize=FS_ANNOT)
        end
        hlines!(ax1, [1.0]; color=:red, linestyle=:dash, linewidth=2)

        # Row 2: polytomy (simulated) null — only if present
        if has_polytomy
            ax2 = Axis(fig[2, j];
                       title="$lbl — polytomy (provisional)",
                       titlesize=FS_TITLE,
                       xlabel="p-value", xlabelsize=FS_LABEL,
                       ylabel="density",  ylabelsize=FS_LABEL,
                       xticklabelsize=FS_TICK, yticklabelsize=FS_TICK,
                       limits=((0, 1), nothing))
            if !isempty(sim_v)
                hist!(ax2, sim_v; normalization=:pdf, bins=20,
                      color=(c, 0.35), strokewidth=1, strokecolor=c)
            else
                text!(ax2, 0.5, 0.5; text="no data", align=(:center, :center),
                      space=:relative, fontsize=FS_ANNOT)
            end
            hlines!(ax2, [1.0]; color=:red, linestyle=:dash, linewidth=2)
        end
    end

    # Legend
    elem_ref = LineElement(color=:red, linestyle=:dash, linewidth=2)
    Legend(fig[nrows + 1, :], [elem_ref], ["Uniform null (density = 1)"];
           orientation=:horizontal, framevisible=false, labelsize=FS_LEGEND)

    out = joinpath(OUT_DIR, "splitting_null_pvalue_distribution.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
