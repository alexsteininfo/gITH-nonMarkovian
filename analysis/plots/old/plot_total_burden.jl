using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Serialization
using CairoMakie
using Statistics

include(joinpath(dirname(@__DIR__), "helpers", "types.jl"))

const MODES = [
    (name = "growth",     label = "Additive (fixed s)",      datadir = "growth",           pfx = "growth_s"),
    (name = "multi",      label = "Multiplicative (fixed)",  datadir = "growth_multi",      pfx = "growth_multi_s"),
    (name = "multi_rand", label = "Multiplicative (random)", datadir = "growth_multi_rand", pfx = "growth_multi_rand_s"),
    (name = "rand",       label = "Random max (Exp)",        datadir = "growth_rand",       pfx = "growth_rand_s"),
]

const S_VALUES  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
const RAW_ROOT  = joinpath(@__DIR__, "..", "..", "data", "raw")
const PLOT_ROOT = joinpath(@__DIR__, "..", "..", "plots", "scMB")

const COL_BURDEN  = (:steelblue, 0.6)
const COL_FITNESS = (:firebrick, 0.6)

const FS_HEAD   = 11
const FS_TITLE  = 10
const FS_LABEL  =  9
const FS_TICK   =  8
const FS_LEGEND =  8

nrows_ea, ncols = 2, 3

function total_unique_mutations(root::BinaryNode{NonMarkovCell})::Int64
    count = root.data.mutations
    if !isnothing(root.left)
        count += total_unique_mutations(root.left)
    end
    if !isnothing(root.right)
        count += total_unique_mutations(root.right)
    end
    return count
end

for mode in MODES
    println("=== Total burden: $(mode.label) ===")
    outdir = joinpath(PLOT_ROOT, mode.name)
    mkpath(outdir)

    fig = Figure(size = (ncols * 260, nrows_ea * 2 * 210 + 100), figure_padding = 12)
    Label(fig[1, 1:ncols]; text = "Total unique driver mutations",
          fontsize = FS_HEAD, font = :bold, tellwidth = false, padding = (0,0,4,0))
    Label(fig[nrows_ea + 2, 1:ncols]; text = "Mean fitness at endpoint",
          fontsize = FS_HEAD, font = :bold, tellwidth = false, padding = (0,8,4,0))

    burden_axes  = Vector{Axis}(undef, length(S_VALUES))
    fitness_axes = Vector{Axis}(undef, length(S_VALUES))

    for (idx, s) in enumerate(S_VALUES)
        r = (idx - 1) ÷ ncols + 1
        c = (idx - 1) % ncols + 1

        ax_b = Axis(fig[r + 1, c];
            title          = "s = $s",
            titlesize      = FS_TITLE,
            xlabel         = r == nrows_ea ? "unique driver mutations" : "",
            ylabel         = c == 1 ? "simulations" : "",
            xlabelsize     = FS_LABEL,
            ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,
            yticklabelsize = FS_TICK,
            spinewidth     = 0.8,
        )
        ax_f = Axis(fig[r + nrows_ea + 2, c];
            title          = "s = $s",
            titlesize      = FS_TITLE,
            xlabel         = r == nrows_ea ? "mean fitness" : "",
            ylabel         = c == 1 ? "simulations" : "",
            xlabelsize     = FS_LABEL,
            ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,
            yticklabelsize = FS_TICK,
            spinewidth     = 0.8,
        )
        burden_axes[idx]  = ax_b
        fitness_axes[idx] = ax_f

        infile = joinpath(RAW_ROOT, mode.datadir, "$(mode.pfx)$(s)_k5.0.jls")
        if !isfile(infile)
            for ax in (ax_b, ax_f)
                text!(ax, 0.5, 0.5; text = "data not found", align = (:center, :center),
                      space = :relative, fontsize = 9, color = :grey60)
            end
            continue
        end

        sims = deserialize(infile)
        println("  s=$s: $(length(sims)) sims")

        burdens   = Int[]
        fitnesses = Float64[]
        burdens   = Int64[total_unique_mutations(something(s.tree_root))
                          for s in sims if !isnothing(s.tree_root)]
        fitnesses = Float64[s.trajectory[end].mean_fitness
                            for s in sims if !isempty(s.trajectory)]

        if !isempty(burdens)
            hist!(ax_b, Float64.(burdens); bins = 30, color = COL_BURDEN)
            text!(ax_b, 0.97, 0.97;
                text     = "mean = $(round(mean(burdens), digits=1))\nvar  = $(round(var(Float64.(burdens)), digits=1))",
                align    = (:right, :top), space = :relative, fontsize = 7.0, color = :black)
        end
        if !isempty(fitnesses)
            hist!(ax_f, fitnesses; bins = 30, color = COL_FITNESS)
            text!(ax_f, 0.97, 0.97;
                text     = "mean = $(round(mean(fitnesses), digits=2))\nvar  = $(round(var(fitnesses), digits=3))",
                align    = (:right, :top), space = :relative, fontsize = 7.0, color = :black)
        end
    end

    linkxaxes!(burden_axes...)
    linkyaxes!(burden_axes...)
    linkxaxes!(fitness_axes...)
    linkyaxes!(fitness_axes...)

    rowgap!(fig.layout, 1,            2.0)
    rowgap!(fig.layout, nrows_ea + 1, 10.0)
    rowgap!(fig.layout, nrows_ea + 2, 2.0)

    outfile = joinpath(outdir, "total_burden_$(mode.name).png")
    save(outfile, fig)
    println("  → $outfile")
end
