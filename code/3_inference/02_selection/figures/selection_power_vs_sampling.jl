# Figure 15: selection_power_vs_sampling.png
#
# Power of the splitting_null test vs. selection coefficient s, stratified by
# sampling depth (n = N / 0.1N / 0.01N). For each timing model, three overlaid
# curves show how detection power degrades as cell-sample size decreases.
# A neutral Type-I error line (from full-tree neutral data) is shown for
# reference.
#
# Inputs:
#   data/inference/selection/selection_1/full_trees/*.jls      (n = N)
#   data/inference/selection/selection_1/subsampled_trees/*.jls  (n < N)
#   data/inference/selection/neutral/full_trees/*.jls           (Type I rate)
# Outputs: figures/3_inference/selection/selection_power_vs_sampling.png
#
# If neither directory is present, the script writes an empty PNG and exits 0.
#
# Usage:
#   julia --project=. code/3_inference/02_selection/figures/selection_power_vs_sampling.jl

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

const SEL1_FULL_DIR = joinpath(DATA, "inference", "selection", "selection_1", "full_trees")
const SEL1_SUB_DIR  = joinpath(DATA, "inference", "selection", "selection_1", "subsampled_trees")
const NEUTRAL_DIR   = joinpath(DATA, "inference", "selection", "neutral", "full_trees")
const OUT_DIR       = joinpath(FIGURES, "3_inference", "selection")
mkpath(OUT_DIR)

const ALPHA = 0.05

# ---------------------------------------------------------------------------
# Stem parsers
# ---------------------------------------------------------------------------

# sel1 full-tree stem: sel1_<timing>_N<N>[_d<d>][_k<k>]_s<s>
function parse_stem_sel1_full(stem::String)
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
    return (; timing, N, d, k, s=sel, n=N)  # n = N for full tree
end

# sel1 subsampled stem: sel1_<timing>_N<N>[_d<d>][_k<k>]_s<s>_n<n>
function parse_stem_sel1_sub(stem::String)
    m = match(r"^(.+)_n(\d+)$", stem)
    isnothing(m) && error("cannot parse subsampled stem: $stem")
    full_stem = m.captures[1]
    n = parse(Int, m.captures[2])
    pf = parse_stem_sel1_full(full_stem)
    return (pf..., n=n)
end

# neutral stem: neutral_<timing>_N<N>[_d<d>[_k<k>]]
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

# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

# Returns Dict: timing -> Vector of (s_val, is_significant, n)
function collect_sel1(dir::String, parse_stem_fn)
    by_timing = Dict{String, Vector{Tuple{Float64, Bool, Int}}}()
    isdir(dir) || return by_timing
    for f in sort(readdir(dir))
        endswith(f, ".jls") || continue
        stem = replace(f, ".jls" => "")
        pf = try parse_stem_fn(stem) catch; continue end
        timing = pf.timing
        n_val  = pf.n
        data = try deserialize(joinpath(dir, f)) catch e; @info "skip $f: $e"; continue end
        rows = get!(by_timing, timing, Tuple{Float64, Bool, Int}[])
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
                push!(rows, (s_val, pv < ALPHA, n_val))
            catch e
                @info "skip sim in $stem: $e"
            end
        end
    end
    return by_timing
end

# Returns Dict: timing -> Float64  (Type I error rate from neutral data)
function collect_neutral_type1(dir::String)
    by_timing = Dict{String, Float64}()
    isdir(dir) || return by_timing
    accum_sig   = Dict{String, Int}()
    accum_total = Dict{String, Int}()
    for f in sort(readdir(dir))
        endswith(f, ".jls") || continue
        stem = replace(f, ".jls" => "")
        pf = try parse_stem_neutral(stem) catch; continue end
        timing = pf.timing
        data = try deserialize(joinpath(dir, f)) catch e; @info "skip $f: $e"; continue end
        for d in data
            isnothing(d) && continue
            try
                tbl = d[:splitting_null_pvalues]
                (isnothing(tbl) || tbl isa Exception) && continue
                pvals  = Float64.(tbl[:pvalue])
                exacts = tbl[:exact]
                for i in eachindex(pvals)
                    Bool(exacts[i]) || continue   # binary (exact) null only
                    isfinite(pvals[i]) || continue
                    accum_total[timing] = get(accum_total, timing, 0) + 1
                    pvals[i] < ALPHA && (accum_sig[timing] = get(accum_sig, timing, 0) + 1)
                end
            catch e
                @info "skip sim in $stem: $e"
            end
        end
    end
    for timing in keys(accum_total)
        tot = accum_total[timing]
        tot > 0 && (by_timing[timing] = get(accum_sig, timing, 0) / tot)
    end
    return by_timing
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    full_data = collect_sel1(SEL1_FULL_DIR, parse_stem_sel1_full)
    sub_data  = collect_sel1(SEL1_SUB_DIR,  parse_stem_sel1_sub)
    type1_data = collect_neutral_type1(NEUTRAL_DIR)

    # Merge by timing
    all_timings = union(keys(full_data), keys(sub_data), keys(type1_data))

    if isempty(all_timings)
        @info "No data found — producing empty figure."
        fig = Figure(size=(900, 400))
        Label(fig[1, 1], "No data (selection_1 dirs empty)"; fontsize=FS_LABEL)
        out = joinpath(OUT_DIR, "selection_power_vs_sampling.png")
        save(out, fig)
        println("wrote (empty) ", out)
        return
    end

    # Combine full and subsampled rows per timing
    merged = Dict{String, Vector{Tuple{Float64, Bool, Int}}}()
    for (timing, rows) in full_data
        merged[timing] = copy(rows)
    end
    for (timing, rows) in sub_data
        lst = get!(merged, timing, Tuple{Float64, Bool, Int}[])
        append!(lst, rows)
    end

    timing_order  = ["deterministic", "gamma", "markov"]
    timing_labels = Dict("deterministic" => "Deterministic",
                         "gamma"         => "Gamma (k=5)",
                         "markov"        => "Markov")
    timing_colors = Dict("deterministic" => COL_DET,
                         "gamma"         => COL_GAMMA,
                         "markov"        => COL_MARKOV)

    timings = filter(t -> haskey(merged, t) || haskey(type1_data, t), timing_order)
    if isempty(timings)
        timings = sort(collect(all_timings))
    end

    fig = Figure(size=(420 * max(1, length(timings)), 560))
    Label(fig[0, :], "Splitting-null power vs. s — by sampling depth  (α = $ALPHA)";
          fontsize=FS_HEAD, halign=:center)

    # Determine the set of n-labels to show consistently across panels.
    # We key n-values as integers; n==N is recorded as pf.N (full stem) so we
    # cannot directly label it without the timing's N.  Instead we annotate via
    # the curve label "full (n=N)" and the numeric n for sub.
    # Assign colours by fraction-of-N bucket: full → panel colour; 0.1N → lighter; 0.01N → lightest
    n_alphas = [(1.0, "full (n=N)"), (0.65, "n = 0.1N"), (0.35, "n = 0.01N")]

    for (col, timing) in enumerate(timings)
        c_base = timing_colors[timing]
        lbl    = get(timing_labels, timing, timing)
        ax = Axis(fig[1, col];
                  title=lbl, titlesize=FS_TITLE,
                  xlabel="selection coefficient s", xlabelsize=FS_LABEL,
                  ylabel="fraction p < $ALPHA",     ylabelsize=FS_LABEL,
                  xticklabelsize=FS_TICK, yticklabelsize=FS_TICK,
                  limits=(nothing, (0.0, 1.05)))

        rows = get(merged, timing, Tuple{Float64, Bool, Int}[])

        if !isempty(rows)
            # All unique n values: classify as "full" (largest n per timing),
            # 0.1N, 0.01N by rank.  Since full-tree n = pf.N = N, the
            # largest n is the full-tree entry.
            all_n = sort(unique(r[3] for r in rows); rev=true)
            n_labels = Dict{Int, String}()
            for (i, nv) in enumerate(all_n)
                if i <= length(n_alphas)
                    n_labels[nv] = n_alphas[i][2]
                else
                    n_labels[nv] = "n=$nv"
                end
            end

            for (i, nv) in enumerate(all_n)
                alpha_val = i <= length(n_alphas) ? n_alphas[i][1] : 0.25
                curve_label = get(n_labels, nv, "n=$nv")
                subset = [(r[1], r[2]) for r in rows if r[3] == nv]
                unique_s = sort(unique(r[1] for r in subset))
                isempty(unique_s) && continue
                powers = Float64[]
                for sv in unique_s
                    hits = [r[2] for r in subset if r[1] == sv]
                    push!(powers, isempty(hits) ? 0.0 : mean(hits))
                end
                c = (c_base, alpha_val)
                lines!(ax, collect(1:length(unique_s)), powers;
                       color=c, linewidth=2, label=curve_label)
                scatter!(ax, collect(1:length(unique_s)), powers;
                         color=c, markersize=8)
                ax.xticks = (collect(1:length(unique_s)),
                             ["s=$sv" for sv in unique_s])
            end
        end

        # Type I error from neutral full trees
        t1 = get(type1_data, timing, NaN)
        if isfinite(t1)
            hlines!(ax, [t1]; color=:gray, linestyle=:dash, linewidth=2,
                    label="Type I (neutral, full)")
        end
        hlines!(ax, [ALPHA]; color=:red, linestyle=:dot, linewidth=1.5,
                label="α = $ALPHA")
    end

    # Legend
    elem_full  = LineElement(color=:black, linewidth=2)
    elem_sub1  = LineElement(color=(:black, 0.65), linewidth=2)
    elem_sub2  = LineElement(color=(:black, 0.35), linewidth=2)
    elem_t1    = LineElement(color=:gray,  linestyle=:dash, linewidth=2)
    elem_alpha = LineElement(color=:red,   linestyle=:dot,  linewidth=1.5)
    Legend(fig[2, :],
           [elem_full, elem_sub1, elem_sub2, elem_t1, elem_alpha],
           ["full (n=N)", "n = 0.1N", "n = 0.01N",
            "Type I error (neutral, full)", "α = $ALPHA"];
           orientation=:horizontal, framevisible=false, labelsize=FS_LEGEND)

    out = joinpath(OUT_DIR, "selection_power_vs_sampling.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
