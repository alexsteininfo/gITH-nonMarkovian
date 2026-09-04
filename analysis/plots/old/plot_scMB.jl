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

const COL_MEAN   = (:teal,      1.0)
const COL_BAND   = (:teal,      0.2)
const COL_SINGLE = (:darkgreen, 1.0)

const FS_HEAD   = 11
const FS_TITLE  = 10
const FS_LABEL  =  9
const FS_TICK   =  8
const FS_LEGEND =  8

nrows_ea, ncols = 2, 3

function mb_histogram(mb::Vector{Int64}, min_j::Int, max_j::Int)
    hist = zeros(Float64, max_j - min_j + 1)
    for j in mb
        min_j <= j <= max_j && (hist[j - min_j + 1] += 1.0)
    end
    return hist
end

for mode in MODES
    println("=== Driver scMB: $(mode.label) ===")
    outdir = joinpath(PLOT_ROOT, mode.name)
    mkpath(outdir)

    mpc_by_s = Vector{Union{Nothing, Vector{Vector{Int64}}}}(undef, length(S_VALUES))
    for (idx, s) in enumerate(S_VALUES)
        f = joinpath(RAW_ROOT, mode.datadir, "$(mode.pfx)$(s)_k5.0.jls")
        if isfile(f)
            sims = deserialize(f)
            mpc_by_s[idx] = [
                isnothing(sim.tree_root) ? Int[] : mutations_per_cell(sim.tree_root)
                for sim in sims
            ]
            println("  s=$s: $(length(sims)) sims")
        else
            mpc_by_s[idx] = nothing
        end
    end

    fig = Figure(size = (ncols * 260, nrows_ea * 220 + 100), figure_padding = 12)
    Label(fig[1, 1:ncols]; text = "Driver mutations per cell — $(mode.label)",
          fontsize = FS_HEAD, font = :bold, tellwidth = false, padding = (0,0,4,0))

    axes = Vector{Axis}(undef, length(S_VALUES))

    for (idx, s) in enumerate(S_VALUES)
        r = (idx - 1) ÷ ncols + 1
        c = (idx - 1) % ncols + 1

        ax = Axis(fig[r + 1, c];
            title          = "s = $s",
            titlesize      = FS_TITLE,
            xlabel         = r == nrows_ea ? "driver mutations per cell, j" : "",
            ylabel         = c == 1 ? "Mⱼ  (cells)" : "",
            xlabelsize     = FS_LABEL,
            ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,
            yticklabelsize = FS_TICK,
            spinewidth     = 0.8,
        )
        axes[idx] = ax

        mpc_list = mpc_by_s[idx]
        if isnothing(mpc_list)
            text!(ax, 0.5, 0.5; text = "data not found", align = (:center, :center),
                  space = :relative, fontsize = 9, color = :grey60)
            continue
        end

        valid = filter(!isempty, mpc_list)
        if isempty(valid)
            text!(ax, 0.5, 0.5; text = "no driver mutations", align = (:center, :center),
                  space = :relative, fontsize = 9, color = :grey60)
            continue
        end

        min_j   = minimum(minimum.(valid))
        max_j   = maximum(maximum.(valid))
        xs      = Float64.(min_j:max_j)
        hists   = [mb_histogram(mb, min_j, max_j) for mb in valid]
        mat     = hcat(hists...)
        mb_mean = vec(mean(mat; dims = 2))
        mb_std  = vec(std(mat;  dims = 2))

        band!(ax,    xs, max.(mb_mean .- mb_std, 0.0), mb_mean .+ mb_std; color = COL_BAND)
        lines!(ax,   xs, mb_mean;  color = COL_MEAN,   linewidth = 2.0)
        scatter!(ax, xs, hists[1]; color = COL_SINGLE, markersize = 3.0)
    end

    linkxaxes!(axes...)
    linkyaxes!(axes...)
    rowgap!(fig.layout, 1, 2.0)

    Legend(fig[nrows_ea + 2, 1:ncols],
        [
            PolyElement(color = COL_BAND),
            LineElement(color = COL_MEAN,   linewidth = 2.0),
            MarkerElement(color = COL_SINGLE, marker = :circle, markersize = 6),
        ],
        ["± std", "mean across simulations", "single simulation"];
        orientation  = :horizontal,
        labelsize    = FS_LEGEND,
        framevisible = false,
        tellwidth    = false,
        patchsize    = (16.0, 8.0),
    )

    outfile = joinpath(outdir, "scMB_driver_$(mode.name).png")
    save(outfile, fig)
    println("  → $outfile")
end
