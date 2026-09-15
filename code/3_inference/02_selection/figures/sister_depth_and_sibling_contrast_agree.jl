# Figure 14: sister_depth_and_sibling_contrast_agree.png
#
# For every driver clade across sel_1 sims: scatter its sibling contrast
# (contrast from selection output) vs. its sister depth ratio (rate_ratio from
# mutation-rate output for the matched stem). Overlay the same scatter for a
# matched sample of neutral clades drawn from neutral selection + mutation-rate
# outputs.
#
# The two measures are independent observables of the same underlying fitness
# signal. Agreement between them validates both.
#
# Inputs:  data/inference/selection/selection_1/full_trees/*.jls  (sibling_contrasts)
#          data/inference/mutationrate/selection_1/full_trees/*.jls (sister_depth_ratios)
#            [NOTE: falls back gracefully if mutation-rate sel_1 directory is absent]
#          data/inference/selection/neutral/full_trees/*.jls        (neutral contrasts)
#          data/inference/mutationrate/neutral/full_trees/*.jls     (neutral ratios)
#
# Outputs: figures/3_inference/selection/sister_depth_and_sibling_contrast_agree.png
#
# Usage:
#   julia --project=. code/3_inference/02_selection/figures/sister_depth_and_sibling_contrast_agree.jl

const ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
using Pkg; Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using EvoTracer
using Serialization, Statistics, Random
using CairoMakie
include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "types_selection1.jl"))
include(joinpath(HELPERS, "types_subsampled.jl"))
include(joinpath(HELPERS, "plotting_functions.jl"))

const SEL1_SEL_DIR  = joinpath(DATA, "inference", "selection", "selection_1", "full_trees")
const SEL1_MUT_DIR  = joinpath(DATA, "inference", "mutationrate", "selection_1", "full_trees")
const NEUT_SEL_DIR  = joinpath(DATA, "inference", "selection", "neutral", "full_trees")
const NEUT_MUT_DIR  = joinpath(DATA, "inference", "mutationrate", "neutral", "full_trees")
const OUT_DIR       = joinpath(FIGURES, "3_inference", "selection")
mkpath(OUT_DIR)

# Load driver (contrast, rate_ratio) pairs by matching sel vs mutation-rate shards.
# Returns Vector of (contrast::Float64, rate_ratio::Float64).
function load_driver_pairs(sel_dir::String, mut_dir::String; is_sel1::Bool)
    pairs = Tuple{Float64, Float64}[]
    isdir(sel_dir) || return pairs
    for f in sort(readdir(sel_dir))
        endswith(f, ".jls") || continue
        stem = replace(f, ".jls" => "")
        sel_data = try deserialize(joinpath(sel_dir, f)) catch e; @info "skip sel $f: $e"; continue end

        # Try to load the matching mutation-rate shard (same stem)
        mut_path = joinpath(mut_dir, "$stem.jls")
        mut_data = if isdir(mut_dir) && isfile(mut_path)
            try deserialize(mut_path) catch e; @info "skip mut $f: $e"; nothing end
        else
            nothing
        end

        for (i, d) in enumerate(sel_data)
            isnothing(d) && continue
            try
                # Get contrast from selection output
                if is_sel1
                    inj = d[:injection]
                    (isnothing(inj) || inj isa Exception) && continue
                    driver_id = Int(inj.driver_cell_id)
                    sc_tbl = d[:sibling_contrasts]
                    (isnothing(sc_tbl) || sc_tbl isa Exception) && continue
                    nodes = Int.(sc_tbl[:node])
                    driver_row_sc = findfirst(==(driver_id), nodes)
                    isnothing(driver_row_sc) && continue
                    contrast = Float64(sc_tbl[:contrast][driver_row_sc])
                    isfinite(contrast) || continue

                    # Get rate_ratio from mutation-rate output — same sim index
                    rate_ratio = NaN
                    if mut_data !== nothing && i <= length(mut_data) && !isnothing(mut_data[i])
                        md = mut_data[i]
                        try
                            mutations = md[:mutations]
                            (isnothing(mutations) || mutations isa Exception) && continue
                            sdr_tbl = mutations[:sister_depth_ratios]
                            (isnothing(sdr_tbl) || sdr_tbl isa Exception) && continue
                            col = :rate_ratio in propertynames(sdr_tbl.data) ? :rate_ratio : :ratio
                            sdr_nodes = Int.(sdr_tbl[:node])
                            driver_row_sdr = findfirst(==(driver_id), sdr_nodes)
                            if !isnothing(driver_row_sdr)
                                rr = Float64(sdr_tbl[col][driver_row_sdr])
                                isfinite(rr) && rr > 0 && (rate_ratio = rr)
                            end
                        catch
                        end
                    end
                    isnan(rate_ratio) && continue
                    push!(pairs, (contrast, rate_ratio))
                end
            catch e
                @info "skip sim in $stem[$i]: $e"
            end
        end
    end
    return pairs
end

# Load neutral (contrast, rate_ratio) pairs: match each sim's any-clade contrast with
# the corresponding rate_ratio from the mutation-rate shard.
function load_neutral_pairs(sel_dir::String, mut_dir::String; max_per_shard::Int = 200)
    pairs = Tuple{Float64, Float64}[]
    isdir(sel_dir) || return pairs
    rng = MersenneTwister(42)
    for f in sort(readdir(sel_dir))
        endswith(f, ".jls") || continue
        stem = replace(f, ".jls" => "")
        sel_data = try deserialize(joinpath(sel_dir, f)) catch e; @info "skip sel $f: $e"; continue end
        mut_path = joinpath(mut_dir, "$stem.jls")
        mut_data = if isdir(mut_dir) && isfile(mut_path)
            try deserialize(mut_path) catch e; nothing end
        else
            nothing
        end

        for (i, d) in enumerate(sel_data)
            isnothing(d) && continue
            try
                sc_tbl = d[:sibling_contrasts]
                (isnothing(sc_tbl) || sc_tbl isa Exception) && continue
                contrasts = Float64.(sc_tbl[:contrast])
                nodes = Int.(sc_tbl[:node])
                good_rows = findall(isfinite, contrasts)
                isempty(good_rows) && continue

                # Rate ratios from mut_data — if unavailable, skip
                mut_d = (mut_data !== nothing && i <= length(mut_data)) ? mut_data[i] : nothing
                isnothing(mut_d) && continue
                mutations = mut_d[:mutations]
                (isnothing(mutations) || mutations isa Exception) && continue
                sdr_tbl = mutations[:sister_depth_ratios]
                (isnothing(sdr_tbl) || sdr_tbl isa Exception) && continue
                col = :rate_ratio in propertynames(sdr_tbl.data) ? :rate_ratio : :ratio
                sdr_nodes = Int.(sdr_tbl[:node])
                sdr_ratios = Float64.(sdr_tbl[col])

                # Build a lookup: node -> rate_ratio
                rr_map = Dict(zip(sdr_nodes, sdr_ratios))

                # Sample min(max_per_shard, good_rows) pairs
                sample_rows = length(good_rows) <= max_per_shard ?
                    good_rows : shuffle(rng, good_rows)[1:max_per_shard]
                for r in sample_rows
                    nd = nodes[r]
                    rr = get(rr_map, nd, NaN)
                    isfinite(rr) && rr > 0 || continue
                    push!(pairs, (contrasts[r], rr))
                end
            catch e
                @info "skip neutral sim in $stem[$i]: $e"
            end
        end
    end
    return pairs
end

function main()
    driver_pairs  = load_driver_pairs(SEL1_SEL_DIR, SEL1_MUT_DIR; is_sel1=true)
    neutral_pairs = load_neutral_pairs(NEUT_SEL_DIR, NEUT_MUT_DIR)

    if isempty(driver_pairs) && isempty(neutral_pairs)
        @info "No data found — producing empty figure."
        fig = Figure(size=(700, 480))
        Label(fig[1, 1], "No data (IN_DIRs empty)"; fontsize=FS_LABEL)
        out = joinpath(OUT_DIR, "sister_depth_and_sibling_contrast_agree.png")
        save(out, fig)
        println("wrote (empty) ", out)
        return
    end

    fig = Figure(size=(700, 580))
    Label(fig[0, 1], "Sister depth ratio vs. sibling contrast: driver vs. neutral";
          fontsize=FS_HEAD, halign=:center)

    ax = Axis(fig[1, 1];
              xlabel="sibling contrast δ_tot", xlabelsize=FS_LABEL,
              ylabel="sister depth rate_ratio (log scale)", ylabelsize=FS_LABEL,
              xticklabelsize=FS_TICK, yticklabelsize=FS_TICK,
              yscale=log10)

    legend_elements = []
    legend_labels   = String[]

    # Neutral clades background
    if !isempty(neutral_pairs)
        nc_x = [p[1] for p in neutral_pairs]
        nc_y = [p[2] for p in neutral_pairs]
        scatter!(ax, nc_x, nc_y;
                 color=(:gray, 0.3), markersize=5,
                 label="neutral clades")
        push!(legend_elements, MarkerElement(color=(:gray, 0.5), marker=:circle, markersize=8))
        push!(legend_labels, "neutral clades")
    end

    # Driver clades foreground
    if !isempty(driver_pairs)
        dr_x = [p[1] for p in driver_pairs]
        dr_y = [p[2] for p in driver_pairs]
        scatter!(ax, dr_x, dr_y;
                 color=(:firebrick, 0.7), markersize=8,
                 label="driver clades (sel_1)")
        push!(legend_elements, MarkerElement(color=(:firebrick, 0.7), marker=:circle, markersize=8))
        push!(legend_labels, "driver clades (sel_1)")
    end

    # Reference lines
    vlines!(ax, [0.0]; color=:black, linestyle=:dash, linewidth=1.5)
    hlines!(ax, [1.0]; color=:black, linestyle=:dot, linewidth=1.5)

    if !isempty(legend_elements)
        Legend(fig[1, 2], legend_elements, legend_labels;
               framevisible=false, labelsize=FS_LEGEND)
    end

    out = joinpath(OUT_DIR, "sister_depth_and_sibling_contrast_agree.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
