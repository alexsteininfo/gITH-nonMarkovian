# Figure 5: mutation_rate_bracket vs. sampling depth.
#
# For each timing model, matches full-tree shards to their subsampled
# counterparts by stem prefix and draws three side-by-side box columns per
# (N, d) group: full, 0.1N, 0.01N — showing lower (blue) and upper (orange)
# endpoints. A red dashed line marks M_TRUE = 2.0.
#
# Inputs:
#   data/inference/mutationrate/neutral/full_trees/*.jls
#   data/inference/mutationrate/neutral/subsampled_trees/*_n<n>.jls
# Outputs: figures/3_inference/mutationrate/bracket_vs_sampling.png
#
# Usage:
#   julia --project=. code/3_inference/01_mutationrate/figures/bracket_vs_sampling.jl

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

function parse_stem_full(stem::String)
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

# Subsampled stem: <full_stem>_n<n>   e.g. neutral_gamma_N10000_d0.5_k5.0_n1000
function parse_stem_sub(stem::String)
    m = match(r"^(.+)_n(\d+)$", stem)
    isnothing(m) && error("cannot parse subsampled stem: $stem")
    full_stem = m.captures[1]
    n = parse(Int, m.captures[2])
    pf = parse_stem_full(full_stem)
    return (pf..., n=n)
end

# Extract lower / upper bracket values from a deserialized shard.
function extract_bracket(data)
    lowers, uppers = Float64[], Float64[]
    for d in data
        isnothing(d) && continue
        try
            tbl = d[:mutations][:mutation_rate_bracket]
            tbl isa Exception && continue
            bd = tbl.data.bound
            vl = tbl.data.value
            li = findfirst(==("lower"), bd)
            ui = findfirst(==("upper"), bd)
            (isnothing(li) || isnothing(ui)) && continue
            v_lo = vl[li]; v_up = vl[ui]
            isfinite(v_lo) && push!(lowers, v_lo)
            isfinite(v_up) && push!(uppers, v_up)
        catch e
            # skip silently
        end
    end
    return lowers, uppers
end

struct ShardEntry
    timing::String
    N::Int
    d::Float64
    k::Float64
    n::Int          # sample size; 0 = full tree
    lowers::Vector{Float64}
    uppers::Vector{Float64}
end

function collect_data()
    entries = ShardEntry[]

    # Full trees
    if isdir(FULL_DIR)
        for f in sort(readdir(FULL_DIR))
            endswith(f, ".jls") || continue
            stem = replace(f, ".jls" => "")
            pf = try parse_stem_full(stem) catch; continue end
            data = try deserialize(joinpath(FULL_DIR, f)) catch e; @info "skip $f: $e"; continue end
            lo, up = extract_bracket(data)
            push!(entries, ShardEntry(pf.timing, pf.N, pf.d, pf.k, 0, lo, up))
        end
    end

    # Subsampled trees
    if isdir(SUB_DIR)
        for f in sort(readdir(SUB_DIR))
            endswith(f, ".jls") || continue
            stem = replace(f, ".jls" => "")
            pf = try parse_stem_sub(stem) catch; continue end
            data = try deserialize(joinpath(SUB_DIR, f)) catch e; @info "skip $f: $e"; continue end
            lo, up = extract_bracket(data)
            push!(entries, ShardEntry(pf.timing, pf.N, pf.d, pf.k, pf.n, lo, up))
        end
    end

    return entries
end

function main()
    entries = collect_data()

    if isempty(entries)
        @info "No data found — producing empty figure."
        fig = Figure(size=(900, 400))
        Label(fig[1,1], "No data"; fontsize=FS_LABEL)
        save(joinpath(OUT_DIR, "bracket_vs_sampling.png"), fig)
        println("wrote (empty) ", joinpath(OUT_DIR, "bracket_vs_sampling.png"))
        return
    end

    timings = sort(unique(e.timing for e in entries))
    fig = Figure(size=(400 * length(timings), 500))
    Label(fig[0, :], "Bracket vs. sampling depth  (m = $M_TRUE)";
          fontsize=FS_HEAD, halign=:center)

    for (col, timing) in enumerate(timings)
        sub = filter(e -> e.timing == timing, entries)
        ax = Axis(fig[1, col];
                  title=timing, titlesize=FS_TITLE,
                  xlabel="(N, d) / sampling", xlabelsize=FS_LABEL,
                  ylabel="m̂", ylabelsize=FS_LABEL,
                  xticklabelsize=FS_TICK - 2, yticklabelsize=FS_TICK)
        hlines!(ax, [M_TRUE]; color=:red, linestyle=:dash, linewidth=2)

        # Group by (N, d): each group gets three x-positions (full, 0.1N, 0.01N)
        groups = sort(unique((e.N, e.d) for e in sub))
        xtick_pos   = Float64[]
        xtick_labels = String[]

        for (gi, (N, d)) in enumerate(groups)
            g_sub = filter(e -> e.N == N && e.d == d, sub)
            # full = n=0; subsampled sorted descending by n
            full_entries = filter(e -> e.n == 0, g_sub)
            samp_entries = sort(filter(e -> e.n > 0, g_sub); by=e -> -e.n)

            all_n_entries = vcat(full_entries, samp_entries)  # up to 3
            n_slots = length(all_n_entries)
            # centre position for this group
            cx = (gi - 1) * (n_slots + 1) + (n_slots + 1) / 2

            for (si, entry) in enumerate(all_n_entries)
                xpos = (gi - 1) * (n_slots + 1) + si
                if !isempty(entry.lowers)
                    boxplot!(ax, fill(xpos - 0.15, length(entry.lowers)), entry.lowers;
                             width=0.25, color=(:steelblue, 0.65), strokewidth=1)
                end
                if !isempty(entry.uppers)
                    boxplot!(ax, fill(xpos + 0.15, length(entry.uppers)), entry.uppers;
                             width=0.25, color=(:darkorange, 0.65), strokewidth=1)
                end
                label = entry.n == 0 ? "full" : "n=$(entry.n)"
                push!(xtick_pos,   xpos)
                push!(xtick_labels, label)
            end

            # Group separator label
            push!(xtick_pos,   cx)
            push!(xtick_labels, "N=$N\nd=$d")
        end

        # Use only the group-centre ticks for clarity
        group_ticks = [(gi - 1) * 4 + 2.5 for gi in 1:length(groups)]
        group_labels = ["N=$N d=$d" for (N, d) in groups]
        ax.xticks = (group_ticks, group_labels)
        ax.xticklabelrotation = π / 6
    end

    # Legend
    elem_lo  = PolyElement(color=(:steelblue,  0.65))
    elem_up  = PolyElement(color=(:darkorange, 0.65))
    elem_ref = LineElement(color=:red, linestyle=:dash)
    Legend(fig[2, :], [elem_lo, elem_up, elem_ref],
           ["lower endpoint", "upper endpoint", "truth m=$M_TRUE"];
           orientation=:horizontal, framevisible=false, labelsize=FS_LEGEND)

    out = joinpath(OUT_DIR, "bracket_vs_sampling.png")
    save(out, fig)
    println("wrote ", out)
end

isinteractive() || main()
