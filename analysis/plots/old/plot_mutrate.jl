using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Serialization
using CairoMakie
using Statistics

include(joinpath(dirname(@__DIR__), "helpers", "types.jl"))
include(joinpath(dirname(@__DIR__), "helpers", "tree_analysis.jl"))

const MODES = [
    (name = "growth",     label = "Additive (fixed s)",      datadir = "growth",           pfx = "growth_s"),
    (name = "multi",      label = "Multiplicative (fixed)",  datadir = "growth_multi",      pfx = "growth_multi_s"),
    (name = "multi_rand", label = "Multiplicative (random)", datadir = "growth_multi_rand", pfx = "growth_multi_rand_s"),
    (name = "rand",       label = "Random max (Exp)",        datadir = "growth_rand",       pfx = "growth_rand_s"),
]

const S_VALUES  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
const RAW_ROOT  = joinpath(@__DIR__, "..", "..", "data", "raw")
const PLOT_ROOT = joinpath(@__DIR__, "..", "..", "plots", "scMB")

const COL_HIST = (:slateblue, 0.7)
const COL_MEAN = :black
const COL_TRUE = :firebrick
const LW_LINE  = 1.5

const FS_HEAD   = 11
const FS_TITLE  = 10
const FS_LABEL  =  9
const FS_TICK   =  8
const FS_LEGEND =  8

nrows_ea, ncols = 2, 3

function inferred_driverrate(mpc::Vector{Int64})
    isempty(mpc) && return NaN
    m = mean(Float64.(mpc))
    iszero(m) && return NaN
    return var(Float64.(mpc)) / m - 1
end

for mode in MODES
    println("=== Inferred ν: $(mode.label) ===")
    outdir = joinpath(PLOT_ROOT, mode.name)
    mkpath(outdir)

    fig = Figure(size = (ncols * 260, nrows_ea * 220 + 100), figure_padding = 12)
    Label(fig[1, 1:ncols]; text = "Inferred driver mutation rate ν — $(mode.label)",
          fontsize = FS_HEAD, font = :bold, tellwidth = false, padding = (0,0,4,0))

    axes   = Vector{Axis}(undef, length(S_VALUES))
    true_nu = NaN

    for (idx, s) in enumerate(S_VALUES)
        r = (idx - 1) ÷ ncols + 1
        c = (idx - 1) % ncols + 1

        ax = Axis(fig[r + 1, c];
            title          = "s = $s",
            titlesize      = FS_TITLE,
            xlabel         = r == nrows_ea ? "inferred ν" : "",
            ylabel         = c == 1 ? "simulations" : "",
            xlabelsize     = FS_LABEL,
            ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,
            yticklabelsize = FS_TICK,
            spinewidth     = 0.8,
        )
        axes[idx] = ax

        infile = joinpath(RAW_ROOT, mode.datadir, "$(mode.pfx)$(s)_k5.0.jls")
        if !isfile(infile)
            text!(ax, 0.5, 0.5; text = "data not found", align = (:center, :center),
                  space = :relative, fontsize = 9, color = :grey60)
            continue
        end

        sims = deserialize(infile)
        println("  s=$s: $(length(sims)) sims")

        isnan(true_nu) && (true_nu = sims[1].params.nu)

        rates = Float64[
            inferred_driverrate(mutations_per_cell(something(sim.tree_root)))
            for sim in sims if !isnothing(sim.tree_root)
        ]
        rates = filter(!isnan, rates)
        isempty(rates) && continue

        hist!(ax, rates; bins = 30, color = COL_HIST)
        vlines!(ax, [mean(rates)]; color = COL_MEAN, linewidth = LW_LINE, linestyle = :dash)
        isnan(true_nu) || vlines!(ax, [true_nu]; color = COL_TRUE, linewidth = LW_LINE, linestyle = :solid)
        text!(ax, 0.97, 0.97;
            text     = "mean = $(round(mean(rates), digits=2))\nvar  = $(round(var(rates), digits=3))",
            align    = (:right, :top),
            space    = :relative,
            fontsize = 7.0,
            color    = :black,
        )
    end

    linkxaxes!(axes...)
    linkyaxes!(axes...)
    rowgap!(fig.layout, 1, 2.0)

    Legend(fig[nrows_ea + 2, 1:ncols],
        [
            PolyElement(color = COL_HIST),
            LineElement(color = COL_MEAN, linewidth = LW_LINE, linestyle = :dash),
            LineElement(color = COL_TRUE, linewidth = LW_LINE, linestyle = :solid),
        ],
        ["distribution of inferred ν", "mean of estimates", "true ν = $true_nu"];
        orientation  = :horizontal,
        labelsize    = FS_LEGEND,
        framevisible = false,
        tellwidth    = false,
        patchsize    = (16.0, 8.0),
    )

    outfile = joinpath(outdir, "mutrate_$(mode.name).png")
    save(outfile, fig)
    println("  → $outfile")
end
