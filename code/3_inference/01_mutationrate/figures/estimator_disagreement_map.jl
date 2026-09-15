# Figure 6: Estimator disagreement heatmap.
#
# For each (timing, k, d, N) cell, computes the per-shard median of
# |m̂_cal - m̂_singleton| / M_TRUE across simulations and displays a heatmap
# per timing model (x = d, y = N (log-scale)) coloured by median relative
# disagreement.
#
# Inputs:  data/inference/mutationrate/neutral/full_trees/*.jls
# Outputs: figures/3_inference/mutationrate/estimator_disagreement_map.png
#
# Usage:
#   julia --project=. code/3_inference/01_mutationrate/figures/estimator_disagreement_map.jl

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

struct DisagRow
    timing::String
    N::Int
    d::Float64
    k::Float64
    med_disagreement::Float64
end

function collect_data()
    rows = DisagRow[]
    isdir(IN_DIR) || return rows
    for f in sort(readdir(IN_DIR))
        endswith(f, ".jls") || continue
        stem = replace(f, ".jls" => "")
        pf = try parse_stem(stem) catch; continue end
        data = try deserialize(joinpath(IN_DIR, f)) catch e; @info "skip $f: $e"; continue end
        diffs = Float64[]
        for d in data
            isnothing(d) && continue
            try
                sv_cal = d[:mutations][:calibrated_tree_length_rate]
                sv_sin = d[:mutations][:singleton_rate]
                (sv_cal isa Exception || isnothing(sv_cal)) && continue
                (sv_sin isa Exception || isnothing(sv_sin)) && continue
                v_cal = sv_cal.value
                v_sin = sv_sin.value
                (isfinite(v_cal) && isfinite(v_sin)) || continue
                push!(diffs, abs(v_cal - v_sin) / M_TRUE)
            catch e
                @info "skip sim in $stem: $e"
            end
        end
        isempty(diffs) && continue
        push!(rows, DisagRow(pf.timing, pf.N, pf.d, pf.k, median(diffs)))
    end
    return rows
end

function main()
    rows = collect_data()

    fig = Figure(size=(400, 400))
    Label(fig[0, :], "Estimator disagreement  |m̂_cal − m̂_sin| / m";
          fontsize=FS_HEAD, halign=:center)

    if isempty(rows)
        @info "No data found in $IN_DIR — producing empty figure."
        ax = Axis(fig[1,1])
        text!(ax, 0.5, 0.5; text="No data", align=(:center, :center),
              space=:relative, fontsize=FS_LABEL)
        save(joinpath(OUT_DIR, "estimator_disagreement_map.png"), fig)
        println("wrote (empty) ", joinpath(OUT_DIR, "estimator_disagreement_map.png"))
        return
    end

    timings = sort(unique(r.timing for r in rows))
    all_Ns  = sort(unique(r.N for r in rows))
    all_ds  = sort(unique(r.d for r in rows))
    vmax    = maximum(r.med_disagreement for r in rows if isfinite(r.med_disagreement))

    n_panels = length(timings)
    fig2 = Figure(size=(max(350, 300 * n_panels), 450))
    Label(fig2[0, :], "Estimator disagreement  |m̂_cal − m̂_sin| / m";
          fontsize=FS_HEAD, halign=:center)

    for (col, timing) in enumerate(timings)
        sub = filter(r -> r.timing == timing, rows)
        ax = Axis(fig2[1, col];
                  title=timing, titlesize=FS_TITLE,
                  xlabel="death rate d", xlabelsize=FS_LABEL,
                  ylabel="population size N", ylabelsize=FS_LABEL,
                  xticklabelsize=FS_TICK, yticklabelsize=FS_TICK,
                  yscale=log10)

        # Build a matrix (d × N)
        d_idx = Dict(d => i for (i,d) in enumerate(all_ds))
        N_idx = Dict(N => i for (i,N) in enumerate(all_Ns))
        mat   = fill(NaN, length(all_ds), length(all_Ns))
        for r in sub
            di = get(d_idx, r.d, nothing)
            ni = get(N_idx, r.N, nothing)
            (isnothing(di) || isnothing(ni)) && continue
            mat[di, ni] = r.med_disagreement
        end

        if all(isnan, mat)
            text!(ax, 0.5, 0.5; text="no data", align=(:center, :center),
                  space=:relative, fontsize=FS_ANNOT)
        else
            hm = heatmap!(ax, all_ds, all_Ns, mat;
                          colormap=:plasma, colorrange=(0, max(vmax, 1e-6)))
            Colorbar(fig2[1, col + n_panels]; limits=(0, max(vmax, 1e-6)),
                     colormap=:plasma, label="median rel. disagreement",
                     labelsize=FS_LABEL)
        end

        ax.xticks = (all_ds, string.(all_ds))
        ax.yticks = (all_Ns, string.(all_Ns))
    end

    out = joinpath(OUT_DIR, "estimator_disagreement_map.png")
    save(out, fig2)
    println("wrote ", out)
end

isinteractive() || main()
