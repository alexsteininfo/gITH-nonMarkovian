# Figure 11: delta_tot_vs_delta_dir_on_driver.png
#
# For each sel_1 simulation, extracts the driver clade's row from direct_contrasts
# and scatters delta_tot (contrast_total) vs. delta_dir (contrast_direct), coloured
# by selection coefficient s. A y=x reference line shows where direct = total.
# Points cluster below the diagonal because the driver clade contains internal
# imbalance from the driver itself.
#
# Inputs:  data/inference/selection/selection_1/full_trees/*.jls
# Outputs: figures/3_inference/selection/delta_tot_vs_delta_dir_on_driver.png
#
# Usage:
#   julia --project=. code/3_inference/02_selection/figures/delta_tot_vs_delta_dir_on_driver.jl

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

function parse_stem_sel1(stem::String)
    s = replace(stem, r"^sel1_" => "")
    parts = split(s, "_")
    timing = parts[1]
    idx = 2
    N = parse(Int, replace(parts[idx], "N" => "")); idx += 1
    d = 0.0; k = Inf; sel = NaN
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
    # Returns Vector of (timing, s_val, delta_tot, delta_dir)
    rows = Tuple{String, Float64, Float64, Float64}[]
    isdir(IN_DIR) || return rows
    for f in sort(readdir(IN_DIR))
        endswith(f, ".jls") || continue
        stem = replace(f, ".jls" => "")
        pf = try parse_stem_sel1(stem) catch; continue end
        timing = pf.timing
        data = try deserialize(joinpath(IN_DIR, f)) catch e; @info "skip $f: $e"; continue end
        for d in data
            isnothing(d) && continue
            try
                inj = d[:injection]
                (isnothing(inj) || inj isa Exception) && continue
                dc_tbl = d[:direct_contrasts]
                (isnothing(dc_tbl) || dc_tbl isa Exception) && continue

                driver_id = Int(inj.driver_cell_id)
                nodes = Int.(dc_tbl[:node])
                driver_row = findfirst(==(driver_id), nodes)
                isnothing(driver_row) && continue

                delta_tot = Float64(dc_tbl[:contrast_total][driver_row])
                delta_dir = Float64(dc_tbl[:contrast_direct][driver_row])
                (isfinite(delta_tot) && isfinite(delta_dir)) || continue

                s_val = try Float64(d[:params].s) catch; pf.s end
                isfinite(s_val) || continue

                push!(rows, (timing, s_val, delta_tot, delta_dir))
            catch e
                @info "skip sim in $stem: $e"
            end
        end
    end
    return rows
end

function main()
    rows = collect_data()
    if isempty(rows)
        @info "No data found in $IN_DIR — producing empty figure."
        fig = Figure(size=(700, 480))
        Label(fig[1, 1], "No data (IN_DIR empty)"; fontsize=FS_LABEL)
        out = joinpath(OUT_DIR, "delta_tot_vs_delta_dir_on_driver.png")
        save(out, fig)
        println("wrote (empty) ", out)
        return
    end

    timing_order  = ["deterministic", "gamma", "markov"]
    timing_colors = Dict("deterministic" => COL_DET, "gamma" => COL_GAMMA, "markov" => COL_MARKOV)
    timing_labels = Dict("deterministic" => "Deterministic", "gamma" => "Gamma (k=5)", "markov" => "Markov")
    timings = filter(t -> any(r -> r[1] == t, rows), timing_order)

    fig = Figure(size=(380 * length(timings), 480))
    Label(fig[0, :], "Driver clade: δ_tot vs. δ_dir, coloured by s";
          fontsize=FS_HEAD, halign=:center)

    # Collect all unique s values to build a consistent color scale
    all_s = sort(unique(r[2] for r in rows))
    s_colormap = :plasma
    s_min, s_max = (isempty(all_s) ? (0.0, 1.0) : (minimum(all_s), maximum(all_s)))
    s_range = s_max > s_min ? s_max - s_min : 1.0

    for (j, timing) in enumerate(timings)
        tr = filter(r -> r[1] == timing, rows)
        isempty(tr) && continue
        lbl = timing_labels[timing]

        ax = Axis(fig[1, j];
                  title=lbl, titlesize=FS_TITLE,
                  xlabel="δ_tot (contrast_total)", xlabelsize=FS_LABEL,
                  ylabel="δ_dir (contrast_direct)", ylabelsize=FS_LABEL,
                  xticklabelsize=FS_TICK, yticklabelsize=FS_TICK,
                  aspect=1)

        tot_vals = [r[3] for r in tr]
        dir_vals = [r[4] for r in tr]
        s_vals   = [r[2] for r in tr]
        colors   = [(s_val - s_min) / s_range for s_val in s_vals]

        scatter!(ax, tot_vals, dir_vals;
                 color=colors, colormap=s_colormap, colorrange=(0, 1),
                 markersize=7, alpha=0.7)

        # y = x reference line
        lims = (min(minimum(tot_vals), minimum(dir_vals)) * 1.05,
                max(maximum(tot_vals), maximum(dir_vals)) * 1.05)
        lines!(ax, collect(lims), collect(lims);
               color=:black, linestyle=:dash, linewidth=2,
               label="δ_dir = δ_tot")
    end

    # Colorbar for s
    if !isempty(rows)
        Colorbar(fig[1, length(timings) + 1];
                 colormap=s_colormap,
                 limits=(s_min, s_max),
                 label="selection coefficient s",
                 labelsize=FS_LABEL,
                 ticklabelsize=FS_TICK)
    end

    # Ref line legend
    elem_ref = LineElement(color=:black, linestyle=:dash, linewidth=2)
    Legend(fig[2, :], [elem_ref], ["δ_dir = δ_tot"];
           orientation=:horizontal, framevisible=false, labelsize=FS_LEGEND)

    out = joinpath(OUT_DIR, "delta_tot_vs_delta_dir_on_driver.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
