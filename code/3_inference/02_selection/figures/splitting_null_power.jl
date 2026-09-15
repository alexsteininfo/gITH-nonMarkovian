# Figure 12: splitting_null_power.png
#
# Power curve: fraction of driver clades with p < 0.05 in splitting_null_pvalues,
# vs. selection coefficient s, for sel_1 data. Overlaid: Type I error rate from
# neutral clades (expected ~0.05 for binary null). One panel per timing model.
#
# Inputs:  data/inference/selection/selection_1/full_trees/*.jls  (driver power)
#          data/inference/selection/neutral/full_trees/*.jls       (Type I rate)
# Outputs: figures/3_inference/selection/splitting_null_power.png
#
# Usage:
#   julia --project=. code/3_inference/02_selection/figures/splitting_null_power.jl

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

const SEL1_DIR    = joinpath(DATA, "inference", "selection", "selection_1", "full_trees")
const NEUTRAL_DIR = joinpath(DATA, "inference", "selection", "neutral", "full_trees")
const OUT_DIR     = joinpath(FIGURES, "3_inference", "selection")
mkpath(OUT_DIR)

const ALPHA = 0.05

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

# For each sel_1 sim: does driver node have p < ALPHA?
function collect_sel1(dir::String)
    # by_timing => Vector of (s_val, driver_significant)
    by_timing = Dict{String, Vector{Tuple{Float64, Bool}}}()
    isdir(dir) || return by_timing
    for f in sort(readdir(dir))
        endswith(f, ".jls") || continue
        stem = replace(f, ".jls" => "")
        pf = try parse_stem_sel1(stem) catch; continue end
        timing = pf.timing
        data = try deserialize(joinpath(dir, f)) catch e; @info "skip $f: $e"; continue end
        rows = get!(by_timing, timing, Tuple{Float64, Bool}[])
        for d in data
            isnothing(d) && continue
            try
                inj = d[:injection]
                (isnothing(inj) || inj isa Exception) && continue
                tbl = d[:splitting_null_pvalues]
                (isnothing(tbl) || tbl isa Exception) && continue

                driver_id = Int(inj.driver_cell_id)
                nodes = Int.(tbl[:node])
                pvals = Float64.(tbl[:pvalue])
                driver_row = findfirst(==(driver_id), nodes)
                isnothing(driver_row) && continue

                pv = pvals[driver_row]
                isfinite(pv) || continue

                s_val = try Float64(d[:params].s) catch; pf.s end
                isfinite(s_val) || continue
                push!(rows, (s_val, pv < ALPHA))
            catch e
                @info "skip sim in $stem: $e"
            end
        end
    end
    return by_timing
end

# Type I error rate: fraction of neutral clades with p < ALPHA
function collect_neutral_type1(dir::String)
    by_timing = Dict{String, Float64}()
    isdir(dir) || return by_timing
    for f in sort(readdir(dir))
        endswith(f, ".jls") || continue
        stem = replace(f, ".jls" => "")
        pf = try parse_stem_neutral(stem) catch; continue end
        timing = pf.timing
        data = try deserialize(joinpath(dir, f)) catch e; @info "skip $f: $e"; continue end
        total = 0; sig = 0
        for d in data
            isnothing(d) && continue
            try
                tbl = d[:splitting_null_pvalues]
                (isnothing(tbl) || tbl isa Exception) && continue
                pvals = Float64.(tbl[:pvalue])
                exacts = tbl[:exact]
                for i in eachindex(pvals)
                    Bool(exacts[i]) || continue   # only binary (exact) null
                    isfinite(pvals[i]) || continue
                    total += 1
                    pvals[i] < ALPHA && (sig += 1)
                end
            catch e
                @info "skip sim in $stem: $e"
            end
        end
        if total > 0
            prev = get(by_timing, timing, 0.0)
            by_timing[timing] = prev + sig  # accumulate then divide
            by_timing["_total_" * timing] = get(by_timing, "_total_" * timing, 0.0) + total
        end
    end
    # Convert sums to fractions
    out = Dict{String, Float64}()
    for timing in ("deterministic", "gamma", "markov")
        sig_key = timing; tot_key = "_total_" * timing
        if haskey(by_timing, sig_key) && haskey(by_timing, tot_key)
            out[timing] = by_timing[sig_key] / by_timing[tot_key]
        end
    end
    return out
end

function main()
    sel1_data    = collect_sel1(SEL1_DIR)
    type1_data   = collect_neutral_type1(NEUTRAL_DIR)

    if isempty(sel1_data) && isempty(type1_data)
        @info "No data found — producing empty figure."
        fig = Figure(size=(900, 400))
        Label(fig[1, 1], "No data (IN_DIRs empty)"; fontsize=FS_LABEL)
        out = joinpath(OUT_DIR, "splitting_null_power.png")
        save(out, fig)
        println("wrote (empty) ", out)
        return
    end

    timing_order  = ["deterministic", "gamma", "markov"]
    timing_colors = Dict("deterministic" => COL_DET, "gamma" => COL_GAMMA, "markov" => COL_MARKOV)
    timing_labels = Dict("deterministic" => "Deterministic", "gamma" => "Gamma (k=5)", "markov" => "Markov")
    timings = filter(t -> haskey(sel1_data, t) || haskey(type1_data, t), timing_order)

    fig = Figure(size=(380 * max(1, length(timings)), 480))
    Label(fig[0, :], "Splitting null power vs. s  (α = $ALPHA)";
          fontsize=FS_HEAD, halign=:center)

    for (j, timing) in enumerate(timings)
        c = timing_colors[timing]
        lbl = timing_labels[timing]
        ax = Axis(fig[1, j];
                  title=lbl, titlesize=FS_TITLE,
                  xlabel="selection coefficient s", xlabelsize=FS_LABEL,
                  ylabel="fraction p < $ALPHA",     ylabelsize=FS_LABEL,
                  xticklabelsize=FS_TICK, yticklabelsize=FS_TICK,
                  limits=(nothing, (0, 1)))

        # Power curve from sel_1
        rows = get(sel1_data, timing, Tuple{Float64,Bool}[])
        if !isempty(rows)
            unique_s = sort(unique(r[1] for r in rows))
            powers = Float64[]
            for sv in unique_s
                subset = [r[2] for r in rows if r[1] == sv]
                push!(powers, isempty(subset) ? 0.0 : mean(subset))
            end
            lines!(ax, collect(1:length(unique_s)), powers;
                   color=c, linewidth=2, label="power")
            scatter!(ax, collect(1:length(unique_s)), powers;
                     color=c, markersize=8)
            ax.xticks = (collect(1:length(unique_s)), ["s=$(sv)" for sv in unique_s])
        end

        # Type I error from neutral
        t1 = get(type1_data, timing, NaN)
        if isfinite(t1)
            hlines!(ax, [t1]; color=:gray, linestyle=:dash, linewidth=2,
                    label="Type I (neutral)")
        end
        hlines!(ax, [ALPHA]; color=:red, linestyle=:dot, linewidth=1.5,
                label="α = $ALPHA")
    end

    # Legend
    elem_power = LineElement(color=:gray, linewidth=2)
    elem_t1    = LineElement(color=:gray, linestyle=:dash, linewidth=2)
    elem_alpha = LineElement(color=:red, linestyle=:dot, linewidth=2)
    Legend(fig[2, :], [elem_power, elem_t1, elem_alpha],
           ["power (sel_1 driver)", "Type I error (neutral)", "α = $ALPHA"];
           orientation=:horizontal, framevisible=false, labelsize=FS_LEGEND)

    out = joinpath(OUT_DIR, "splitting_null_power.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
