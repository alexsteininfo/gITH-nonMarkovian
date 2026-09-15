# Figure 10: driver_direct_contrast_rank_vs_s.png
#
# For each sel_1 simulation, identifies the driver clade by matching
# injection.driver_cell_id against the direct_contrasts node column,
# then computes the driver's rank (from the top) in the direct contrast
# distribution for that tree. Boxplot of rank vs. selection coefficient s.
# One panel per timing model.
#
# Inputs:  data/inference/selection/selection_1/full_trees/*.jls
# Outputs: figures/3_inference/selection/driver_direct_contrast_rank_vs_s.png
#
# Usage:
#   julia --project=. code/3_inference/02_selection/figures/driver_direct_contrast_rank_vs_s.jl

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

const IN_DIR  = joinpath(DATA, "inference", "selection", "selection_1", "full_trees")
const OUT_DIR = joinpath(FIGURES, "3_inference", "selection")
mkpath(OUT_DIR)

# Parse a sel_1 stem: sel1_<timing>_N<N>_d<d>_k<k>_s<s>
# Deterministic stems omit _d and _k fields.
function parse_stem_sel1(stem::String)
    s = replace(stem, r"^sel1_" => "")
    parts = split(s, "_")
    timing = parts[1]
    idx = 2
    N = parse(Int, replace(parts[idx], "N" => "")); idx += 1
    d = 0.0
    k = Inf
    sel = NaN
    if idx <= length(parts) && startswith(parts[idx], "d")
        d = parse(Float64, replace(parts[idx], "d" => "")); idx += 1
    end
    if idx <= length(parts) && startswith(parts[idx], "k")
        k = parse(Float64, replace(parts[idx], "k" => "")); idx += 1
    end
    if idx <= length(parts) && startswith(parts[idx], "s")
        sel = parse(Float64, replace(parts[idx], "s" => "")); idx += 1
    end
    return (; timing, N, d, k, s=sel)
end

function collect_data()
    # by_timing => Vector of (s_value, rank_fraction)
    by_timing = Dict{String, Vector{Tuple{Float64, Float64}}}()
    isdir(IN_DIR) || return by_timing
    for f in sort(readdir(IN_DIR))
        endswith(f, ".jls") || continue
        stem = replace(f, ".jls" => "")
        pf = try parse_stem_sel1(stem) catch; continue end
        timing = pf.timing
        data = try deserialize(joinpath(IN_DIR, f)) catch e; @info "skip $f: $e"; continue end
        rows = get!(by_timing, timing, Tuple{Float64,Float64}[])
        for d in data
            isnothing(d) && continue
            try
                inj = d[:injection]
                (isnothing(inj) || inj isa Exception) && continue
                dc_tbl = d[:direct_contrasts]
                (isnothing(dc_tbl) || dc_tbl isa Exception) && continue

                driver_id = Int(inj.driver_cell_id)
                nodes = Int.(dc_tbl[:node])
                dir_vals = Float64.(dc_tbl[:contrast_direct])

                # Find the driver's row: match by node id
                driver_row = findfirst(==(driver_id), nodes)
                isnothing(driver_row) && continue

                driver_val = dir_vals[driver_row]
                isfinite(driver_val) || continue

                # Rank from top (rank 1 = best): how many are >= driver_val
                finite_vals = filter(isfinite, dir_vals)
                isempty(finite_vals) && continue
                rank = count(>=(driver_val), finite_vals)  # rank from top (1 = highest)
                rank_frac = rank / length(finite_vals)

                # s value from params or filename
                s_val = try Float64(d[:params].s) catch; pf.s end
                isfinite(s_val) || continue

                push!(rows, (s_val, rank_frac))
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
        out = joinpath(OUT_DIR, "driver_direct_contrast_rank_vs_s.png")
        save(out, fig)
        println("wrote (empty) ", out)
        return
    end

    timing_order  = ["deterministic", "gamma", "markov"]
    timing_colors = Dict("deterministic" => COL_DET, "gamma" => COL_GAMMA, "markov" => COL_MARKOV)
    timing_labels = Dict("deterministic" => "Deterministic", "gamma" => "Gamma (k=5)", "markov" => "Markov")
    timings = filter(t -> haskey(by_timing, t), timing_order)

    fig = Figure(size=(380 * length(timings), 480))
    Label(fig[0, :], "Driver direct contrast rank vs. s (rank 1 = highest δ_dir)";
          fontsize=FS_HEAD, halign=:center)

    for (j, timing) in enumerate(timings)
        rows = by_timing[timing]
        isempty(rows) && continue
        c = timing_colors[timing]
        lbl = timing_labels[timing]

        ax = Axis(fig[1, j];
                  title=lbl, titlesize=FS_TITLE,
                  xlabel="selection coefficient s", xlabelsize=FS_LABEL,
                  ylabel="driver rank (fraction from top)", ylabelsize=FS_LABEL,
                  xticklabelsize=FS_TICK, yticklabelsize=FS_TICK)

        s_vals = [r[1] for r in rows]
        rank_vals = [r[2] for r in rows]
        unique_s = sort(unique(s_vals))

        for (xi, sv) in enumerate(unique_s)
            ys = [rank_vals[i] for i in eachindex(s_vals) if s_vals[i] == sv]
            isempty(ys) && continue
            boxplot!(ax, fill(xi, length(ys)), ys;
                     color=(c, 0.6), strokewidth=1,
                     width=0.6)
        end
        ax.xticks = (collect(1:length(unique_s)), ["s=$(sv)" for sv in unique_s])
        hlines!(ax, [0.1]; color=:red, linestyle=:dash, linewidth=1.5,
                label="top decile")
    end

    # Legend
    elem_ref = LineElement(color=:red, linestyle=:dash, linewidth=2)
    Legend(fig[2, :], [elem_ref], ["top decile (rank fraction = 0.1)"];
           orientation=:horizontal, framevisible=false, labelsize=FS_LEGEND)

    out = joinpath(OUT_DIR, "driver_direct_contrast_rank_vs_s.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
