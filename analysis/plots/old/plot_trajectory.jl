using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Serialization
using CairoMakie

include(joinpath(dirname(@__DIR__), "helpers", "types.jl"))

const MODES = [
    (name = "growth",     label = "Additive (fixed s)",      datadir = "growth",           pfx = "growth_s"),
    (name = "multi",      label = "Multiplicative (fixed)",  datadir = "growth_multi",      pfx = "growth_multi_s"),
    (name = "multi_rand", label = "Multiplicative (random)", datadir = "growth_multi_rand", pfx = "growth_multi_rand_s"),
    (name = "rand",       label = "Random max (Exp)",        datadir = "growth_rand",       pfx = "growth_rand_s"),
]

const S_VALUES     = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
const N_TRAJ       = 1
const DATADIR_ROOT = joinpath(@__DIR__, "..", "..", "data", "raw")
const PLOTDIR      = joinpath(@__DIR__, "..", "..", "plots", "trajectories")

const COL_N    = :black
const COL_F    = :firebrick
const LW       = 2.0
const FS_HEAD  = 11
const FS_TITLE = 10
const FS_LABEL =  9
const FS_TICK  =  8
const FS_LEGEND =  8

mkpath(PLOTDIR)
nrows_ea, ncols = 2, 3

for mode in MODES
    println("=== Trajectory: $(mode.label) ===")

    fig = Figure(size = (ncols * 260, nrows_ea * 2 * 220 + 100), figure_padding = 12)
    Label(fig[1, 1:ncols]; text = "Population size (log scale)",
          fontsize = FS_HEAD, font = :bold, tellwidth = false, padding = (0,0,4,0))
    Label(fig[nrows_ea + 2, 1:ncols]; text = "Mean fitness",
          fontsize = FS_HEAD, font = :bold, tellwidth = false, padding = (0,8,4,0))

    axes_n = Vector{Axis}(undef, length(S_VALUES))
    axes_f = Vector{Axis}(undef, length(S_VALUES))

    for (idx, s) in enumerate(S_VALUES)
        r = (idx - 1) ÷ ncols + 1
        c = (idx - 1) % ncols + 1

        ax_n = Axis(fig[r + 1, c];
            yscale         = log10,
            title          = "s = $s",
            titlesize      = FS_TITLE,
            xlabel         = r == nrows_ea ? "time" : "",
            ylabel         = c == 1 ? "cells" : "",
            xlabelsize     = FS_LABEL,
            ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,
            yticklabelsize = FS_TICK,
            spinewidth     = 0.8,
        )
        ax_f = Axis(fig[r + nrows_ea + 2, c];
            title          = "s = $s",
            titlesize      = FS_TITLE,
            xlabel         = r == nrows_ea ? "time" : "",
            ylabel         = c == 1 ? "mean fitness" : "",
            xlabelsize     = FS_LABEL,
            ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,
            yticklabelsize = FS_TICK,
            spinewidth     = 0.8,
        )
        axes_n[idx] = ax_n
        axes_f[idx] = ax_f

        infile = joinpath(DATADIR_ROOT, mode.datadir, "$(mode.pfx)$(s)_k5.0.jls")
        if !isfile(infile)
            for ax in (ax_n, ax_f)
                text!(ax, 0.5, 0.5; text = "data not found", align = (:center, :center),
                      space = :relative, fontsize = 9, color = :grey60)
            end
            continue
        end

        sims = deserialize(infile)
        println("  s=$s: $(min(N_TRAJ, length(sims))) trajectories")
        for sim in sims[1:min(N_TRAJ, length(sims))]
            traj = sim.trajectory
            ts   = [p.t            for p in traj]
            Ns   = Float64[p.N_total      for p in traj]
            fs   = Float64[p.mean_fitness for p in traj]
            lines!(ax_n, ts, Ns; color = COL_N, linewidth = LW)
            lines!(ax_f, ts, fs; color = COL_F, linewidth = LW)
        end
    end

    linkyaxes!(axes_n...)
    linkyaxes!(axes_f...)
    linkxaxes!(axes_n...)
    linkxaxes!(axes_f...)

    rowgap!(fig.layout, 1,            2.0)
    rowgap!(fig.layout, nrows_ea + 1, 10.0)
    rowgap!(fig.layout, nrows_ea + 2, 2.0)

    Legend(fig[nrows_ea * 2 + 3, 1:ncols],
        [
            LineElement(color = COL_N, linewidth = LW),
            LineElement(color = COL_F, linewidth = LW),
        ],
        ["population size (log₁₀)", "mean fitness"];
        orientation  = :horizontal,
        labelsize    = FS_LEGEND,
        framevisible = false,
        tellwidth    = false,
        patchsize    = (16.0, 8.0),
    )

    outfile = joinpath(PLOTDIR, "trajectory_$(mode.name).png")
    save(outfile, fig)
    println("  → $outfile")
end
