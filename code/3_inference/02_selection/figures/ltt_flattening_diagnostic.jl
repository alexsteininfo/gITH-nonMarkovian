# Figure 13: ltt_flattening_diagnostic.png
#
# Lineages-through-time (LTT) curves for representative neutral and sel_1 shards.
# All simulations from a shard are overlaid with alpha=0.1. A shaded region marks
# the "informative growth phase" where the LTT derivative exceeds 10% of the
# initial slope (i.e. before flattening becomes severe).
#
# Inputs:  data/inference/selection/neutral/full_trees/*.jls
#          data/inference/selection/selection_1/full_trees/*.jls
# Outputs: figures/3_inference/selection/ltt_flattening_diagnostic.png
#
# Usage:
#   julia --project=. code/3_inference/02_selection/figures/ltt_flattening_diagnostic.jl

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
using Pkg; Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using EvoTracer
using Serialization, Statistics
using CairoMakie
include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "plotting_functions.jl"))

const NEUTRAL_DIR = joinpath(DATA, "inference", "selection", "neutral", "full_trees")
const SEL1_DIR    = joinpath(DATA, "inference", "selection", "selection_1", "full_trees")
const OUT_DIR     = joinpath(FIGURES, "3_inference", "selection")
mkpath(OUT_DIR)

const DERIVATIVE_THRESHOLD = 0.10   # fraction of initial slope

# Load LTT curves from one shard file; returns Vector of (dist, count) pairs
function load_ltt_curves(path::String)
    data = try deserialize(path) catch e; @info "skip $(basename(path)): $e"; return [] end
    curves = Tuple{Vector{Float64}, Vector{Int}}[]
    for d in data
        isnothing(d) && continue
        try
            tbl = d[:lineages_through_time]
            (isnothing(tbl) || tbl isa Exception) && continue
            dists  = Float64.(tbl[:distance])
            counts = Int.(tbl[:lineages])
            (isempty(dists) || isempty(counts)) && continue
            push!(curves, (dists, counts))
        catch e
            @info "ltt parse: $e"
        end
    end
    return curves
end

# Compute the inflection distance where the derivative drops below threshold of
# the initial slope, using finite differences on the pooled LTT.
function informative_cutoff(dists::Vector{Float64}, counts::Vector{Int};
                             threshold::Float64 = DERIVATIVE_THRESHOLD)
    length(dists) < 2 && return last(dists)
    diffs  = diff(Float64.(counts))
    gaps   = diff(dists)
    slopes = [gaps[i] > 0 ? diffs[i] / gaps[i] : 0.0 for i in eachindex(gaps)]
    init   = maximum(slopes)
    init <= 0 && return last(dists)
    cutoff_idx = findlast(s -> s >= threshold * init, slopes)
    isnothing(cutoff_idx) && return first(dists)
    return dists[cutoff_idx + 1]
end

function pick_shard(dir::String, timing::String, prefix::String)
    isdir(dir) || return nothing
    for f in sort(readdir(dir))
        endswith(f, ".jls") || continue
        startswith(f, "$(prefix)_$(timing)_") && return joinpath(dir, f)
    end
    return nothing
end

function main()
    timing_order  = ["deterministic", "gamma", "markov"]
    timing_colors = Dict("deterministic" => COL_DET, "gamma" => COL_GAMMA, "markov" => COL_MARKOV)
    timing_labels = Dict("deterministic" => "Deterministic", "gamma" => "Gamma (k=5)", "markov" => "Markov")

    # Gather representative shards: one neutral + one sel_1 per timing model
    panel_specs = Tuple{String, String, String}[]  # (label, path, timing)
    for timing in timing_order
        n_path = pick_shard(NEUTRAL_DIR, timing, "neutral")
        s_path = pick_shard(SEL1_DIR,    timing, "sel1")
        n_path === nothing || push!(panel_specs, ("$timing / neutral", n_path, timing))
        s_path === nothing || push!(panel_specs, ("$timing / sel_1", s_path, timing))
    end

    if isempty(panel_specs)
        @info "No data found — producing empty figure."
        fig = Figure(size=(900, 400))
        Label(fig[1, 1], "No data (IN_DIRs empty)"; fontsize=FS_LABEL)
        out = joinpath(OUT_DIR, "ltt_flattening_diagnostic.png")
        save(out, fig)
        println("wrote (empty) ", out)
        return
    end

    ncols = min(3, length(panel_specs))
    nrows = cld(length(panel_specs), ncols)
    fig = Figure(size=(380 * ncols, 350 * nrows + 60))
    Label(fig[0, :], "LTT diagnostic — shaded = informative growth phase";
          fontsize=FS_HEAD, halign=:center)

    for (idx, (lbl, path, timing)) in enumerate(panel_specs)
        row = cld(idx, ncols)
        col = mod1(idx, ncols)
        c = timing_colors[timing]

        ax = Axis(fig[row, col];
                  title=lbl, titlesize=FS_TITLE,
                  xlabel="molecular distance", xlabelsize=FS_LABEL,
                  ylabel="lineages",           ylabelsize=FS_LABEL,
                  xticklabelsize=FS_TICK, yticklabelsize=FS_TICK)

        curves = load_ltt_curves(path)
        if isempty(curves)
            text!(ax, 0.5, 0.5; text="no data", align=(:center, :center),
                  space=:relative, fontsize=FS_ANNOT)
            continue
        end

        # Plot all curves at alpha=0.1
        for (dists, counts) in curves
            lines!(ax, dists, Float64.(counts);
                   color=(c, 0.1), linewidth=1)
        end

        # Pool the first curve for computing the informative cutoff
        ref_d, ref_c = curves[1]
        cutoff = informative_cutoff(ref_d, ref_c; threshold=DERIVATIVE_THRESHOLD)

        # Shade the informative region using vspan
        x_min = first(ref_d)
        vspan!(ax, x_min, cutoff; color=(:gold, 0.25))
        vlines!(ax, [cutoff]; color=:black, linestyle=:dash, linewidth=1.5,
                label="$(round(Int, DERIVATIVE_THRESHOLD*100))% threshold")
    end

    out = joinpath(OUT_DIR, "ltt_flattening_diagnostic.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
