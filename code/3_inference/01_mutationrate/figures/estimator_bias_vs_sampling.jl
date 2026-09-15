# Figure 16: estimator_bias_vs_sampling.png
#
# Relative bias (m̂ − M_TRUE) / M_TRUE of three mutation-rate estimators vs.
# log10(n) (sample size), stratified by timing model. For each (timing, estimator,
# n) the per-shard median bias is computed; then one line per estimator is drawn
# across n values.
#
# Estimators:
#   mutation_rate_bracket.lower
#   mutation_rate_bracket.upper
#   calibrated_tree_length_rate.value
#
# n = N (full tree) appears at the right of each x-axis; subsampled at 0.1N and
# 0.01N (if present) appear at smaller log10(n) values.
#
# Inputs:
#   data/inference/mutationrate/neutral/full_trees/*.jls
#   data/inference/mutationrate/neutral/subsampled_trees/*_n<n>.jls
# Outputs: figures/3_inference/mutationrate/estimator_bias_vs_sampling.png
#
# If no data are present, the script writes an empty PNG and exits 0.
#
# Usage:
#   julia --project=. code/3_inference/01_mutationrate/figures/estimator_bias_vs_sampling.jl

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

const FULL_DIR = joinpath(DATA, "inference", "mutationrate", "neutral", "full_trees")
const SUB_DIR  = joinpath(DATA, "inference", "mutationrate", "neutral", "subsampled_trees")
const OUT_DIR  = joinpath(FIGURES, "3_inference", "mutationrate")
mkpath(OUT_DIR)

const M_TRUE = 2.0

# Estimator specs: (label_for_plot, extractor_function)
# Each extractor receives a single Dict `d` (one sim) and returns a Float64 or NaN.
const ESTIMATORS = [
    ("bracket.lower",     d -> begin
        tbl = d[:mutations][:mutation_rate_bracket]
        tbl isa Exception && return NaN
        bd = tbl.data.bound; vl = tbl.data.value
        li = findfirst(==("lower"), bd)
        isnothing(li) ? NaN : Float64(vl[li])
    end),
    ("bracket.upper",     d -> begin
        tbl = d[:mutations][:mutation_rate_bracket]
        tbl isa Exception && return NaN
        bd = tbl.data.bound; vl = tbl.data.value
        ui = findfirst(==("upper"), bd)
        isnothing(ui) ? NaN : Float64(vl[ui])
    end),
    ("calib_tree_rate",   d -> begin
        sv = d[:mutations][:calibrated_tree_length_rate]
        (isnothing(sv) || sv isa Exception) && return NaN
        Float64(sv.value)
    end),
]

# ---------------------------------------------------------------------------
# Stem parsers
# ---------------------------------------------------------------------------

function parse_stem_full(stem::String)
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

# Subsampled stem: <full_stem>_n<n>
function parse_stem_sub(stem::String)
    m = match(r"^(.+)_n(\d+)$", stem)
    isnothing(m) && error("cannot parse subsampled stem: $stem")
    full_stem = m.captures[1]
    n = parse(Int, m.captures[2])
    pf = parse_stem_full(full_stem)
    return (pf..., n=n)
end

# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

# Returns Vector of NamedTuples:
#   (timing, N, d, k, n, est_label, med_bias)
# where n = N for full-tree entries, and med_bias = median((m̂ − M_TRUE)/M_TRUE)
# over all sims in the shard.
function collect_data()
    entries = Any[]

    # Full trees
    if isdir(FULL_DIR)
        for f in sort(readdir(FULL_DIR))
            endswith(f, ".jls") || continue
            stem = replace(f, ".jls" => "")
            pf = try parse_stem_full(stem) catch; continue end
            data = try deserialize(joinpath(FULL_DIR, f)) catch e; @info "skip $f: $e"; continue end
            for (elabel, efn) in ESTIMATORS
                vals = Float64[]
                for d in data
                    isnothing(d) && continue
                    v = try efn(d) catch; NaN end
                    isfinite(v) && push!(vals, (v - M_TRUE) / M_TRUE)
                end
                isempty(vals) && continue
                push!(entries, (timing=pf.timing, N=pf.N, d=pf.d, k=pf.k,
                                n=pf.N, est=elabel, med_bias=median(vals)))
            end
        end
    end

    # Subsampled trees
    if isdir(SUB_DIR)
        for f in sort(readdir(SUB_DIR))
            endswith(f, ".jls") || continue
            stem = replace(f, ".jls" => "")
            pf = try parse_stem_sub(stem) catch; continue end
            data = try deserialize(joinpath(SUB_DIR, f)) catch e; @info "skip $f: $e"; continue end
            for (elabel, efn) in ESTIMATORS
                vals = Float64[]
                for d in data
                    isnothing(d) && continue
                    v = try efn(d) catch; NaN end
                    isfinite(v) && push!(vals, (v - M_TRUE) / M_TRUE)
                end
                isempty(vals) && continue
                push!(entries, (timing=pf.timing, N=pf.N, d=pf.d, k=pf.k,
                                n=pf.n, est=elabel, med_bias=median(vals)))
            end
        end
    end

    return entries
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    entries = collect_data()

    if isempty(entries)
        @info "No data found — producing empty figure."
        fig = Figure(size=(900, 400))
        Label(fig[1, 1], "No data (FULL_DIR and SUB_DIR empty)"; fontsize=FS_LABEL)
        out = joinpath(OUT_DIR, "estimator_bias_vs_sampling.png")
        save(out, fig)
        println("wrote (empty) ", out)
        return
    end

    timing_order  = ["deterministic", "gamma", "markov"]
    timing_labels = Dict("deterministic" => "Deterministic",
                         "gamma"         => "Gamma (k=5)",
                         "markov"        => "Markov")
    timing_colors = Dict("deterministic" => COL_DET,
                         "gamma"         => COL_GAMMA,
                         "markov"        => COL_MARKOV)

    timings = filter(t -> any(e -> e.timing == t, entries), timing_order)
    if isempty(timings)
        timings = sort(unique(e.timing for e in entries))
    end

    # Estimator visual encoding: line style
    est_labels  = [e[1] for e in ESTIMATORS]
    est_styles  = Dict(est_labels[1] => :dash,
                       est_labels[2] => :dot,
                       est_labels[3] => :solid)
    est_markers = Dict(est_labels[1] => :circle,
                       est_labels[2] => :rect,
                       est_labels[3] => :diamond)
    est_display = Dict("bracket.lower"   => "bracket lower",
                       "bracket.upper"   => "bracket upper",
                       "calib_tree_rate" => "calib. tree-length rate")

    fig = Figure(size=(420 * max(1, length(timings)), 560))
    Label(fig[0, :], "Estimator relative bias vs. log₁₀(n)  (M_TRUE = $M_TRUE)";
          fontsize=FS_HEAD, halign=:center)

    for (col, timing) in enumerate(timings)
        c_base = timing_colors[timing]
        lbl    = get(timing_labels, timing, timing)
        ax = Axis(fig[1, col];
                  title=lbl, titlesize=FS_TITLE,
                  xlabel="log₁₀(n)", xlabelsize=FS_LABEL,
                  ylabel="median (m̂ − M) / M", ylabelsize=FS_LABEL,
                  xticklabelsize=FS_TICK, yticklabelsize=FS_TICK)

        # Zero-bias reference
        hlines!(ax, [0.0]; color=:black, linestyle=:dash, linewidth=1.5,
                label="zero bias")

        sub = filter(e -> e.timing == timing, entries)
        isempty(sub) && continue

        for (elabel, _) in ESTIMATORS
            esub = filter(e -> e.est == elabel, sub)
            isempty(esub) && continue

            # Aggregate by n: within a timing panel there may be multiple (N, d, k)
            # parameter sets; average their median biases at each n.
            unique_n = sort(unique(e.n for e in esub))
            log_n    = log10.(unique_n)
            agg_bias = Float64[]
            for nv in unique_n
                biases = [e.med_bias for e in esub if e.n == nv]
                push!(agg_bias, mean(biases))
            end

            ls = get(est_styles,  elabel, :solid)
            mk = get(est_markers, elabel, :circle)
            dl = get(est_display, elabel, elabel)
            lines!(ax, log_n, agg_bias;
                   color=c_base, linestyle=ls, linewidth=2, label=dl)
            scatter!(ax, log_n, agg_bias;
                     color=c_base, marker=mk, markersize=9, strokewidth=1,
                     label=dl)

            # Tick labels: show n on x-axis
            ax.xticks = (log_n, ["n=$nv" for nv in unique_n])
            ax.xticklabelrotation = π / 6
        end
    end

    # Legend: one entry per estimator
    elem_lo   = [LineElement(color=:black, linestyle=:dash, linewidth=2),
                 MarkerElement(color=:black, marker=:circle, markersize=9)]
    elem_up   = [LineElement(color=:black, linestyle=:dot, linewidth=2),
                 MarkerElement(color=:black, marker=:rect,   markersize=9)]
    elem_cal  = [LineElement(color=:black, linestyle=:solid, linewidth=2),
                 MarkerElement(color=:black, marker=:diamond, markersize=9)]
    Legend(fig[2, :],
           [elem_lo, elem_up, elem_cal],
           ["bracket lower", "bracket upper", "calib. tree-length rate"];
           orientation=:horizontal, framevisible=false, labelsize=FS_LEGEND)

    out = joinpath(OUT_DIR, "estimator_bias_vs_sampling.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
