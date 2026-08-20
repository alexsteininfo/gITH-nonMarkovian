using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Distributions
using Serialization
using CairoMakie
using Statistics
using Random

include(joinpath(dirname(dirname(@__DIR__)), "helpers", "types.jl"))
include(joinpath(dirname(@__DIR__), "plotting_functions.jl"))

# ── Parameters ────────────────────────────────────────────────────────────────

const b        = 1.0
const k        = 5.0
const D_VALUES = [0.0, 0.5, 0.9]
const N_LOAD   = 10_000   # use largest trees for richest empirical distribution

const RAW = joinpath(@__DIR__, "..", "..", "..", "data", "raw",  "neutral")
const OUT  = joinpath(@__DIR__, "..", "..", "..", "plots", "neutral")

xs = range(0.0, 6.0; length = 400)

# ── Figure factory ────────────────────────────────────────────────────────────
#
# Produces a two-row timing figure:
#   Top row    — theoretical PDFs for birth (and death for d > 0)
#   Bottom row — empirical inter-division times from simulation trees
#
# dist_birth(d)  : returns the birth distribution for a given d
# dist_death(d)  : returns the death distribution (only called for d > 0)
# raw_path(d)    : returns the raw .jls file path for a given d
# label_div      : legend label for the division distribution line
# label_death    : legend label for the death distribution line
# outpath        : where to save the figure

function make_timing_figure(;
    dist_birth,
    dist_death,
    raw_path,
    label_div,
    label_death,
    outpath,
)
    ncols = length(D_VALUES)
    fig   = Figure(size = (ncols * 300, 2 * 240 + 80), figure_padding = 14)

    Label(fig[1, 1:ncols];
        text      = "Implemented (theoretical PDFs)",
        fontsize  = FS_HEAD, font = :bold,
        tellwidth = false, padding = (0, 0, 4, 0))
    Label(fig[3, 1:ncols];
        text      = "Measured (inter-division times from simulation trees)",
        fontsize  = FS_HEAD, font = :bold,
        tellwidth = false, padding = (0, 8, 4, 0))

    axes_theory  = Vector{Axis}(undef, ncols)
    axes_empiric = Vector{Axis}(undef, ncols)

    for (ci, d) in enumerate(D_VALUES)
        bd = dist_birth(d)

        # Row 1: theoretical ─────────────────────────────────────────────────
        ax_t = Axis(fig[2, ci];
            title          = "d = $d",
            titlesize      = FS_TITLE,
            xlabel         = "waiting time",
            ylabel         = ci == 1 ? "probability density" : "",
            xlabelsize     = FS_LABEL, ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,  yticklabelsize = FS_TICK,
            spinewidth     = 0.8)
        axes_theory[ci] = ax_t

        lines!(ax_t, collect(xs), pdf.(bd, xs); color = COL_DIV, linewidth = 2.0)

        if d > 0.0
            dd = dist_death(d)
            lines!(ax_t, collect(xs), pdf.(dd, xs);
                   color = COL_DEATH, linewidth = 2.0, linestyle = :dash)
        end

        # Row 2: empirical ───────────────────────────────────────────────────
        ax_e = Axis(fig[4, ci];
            title          = "d = $d",
            titlesize      = FS_TITLE,
            xlabel         = "inter-division time",
            ylabel         = ci == 1 ? "probability density" : "",
            xlabelsize     = FS_LABEL, ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,  yticklabelsize = FS_TICK,
            spinewidth     = 0.8)
        axes_empiric[ci] = ax_e

        fpath = raw_path(d)
        if !isfile(fpath)
            text!(ax_e, 0.5, 0.5; text = "data not found",
                  align = (:center, :center), space = :relative,
                  fontsize = 9, color = :grey60)
            continue
        end

        sims      = deserialize(fpath)
        lifetimes = vcat([celllifetimes(something(s.tree_root))
                          for s in sims if !isnothing(s.tree_root)]...)
        println("  d=$d: $(length(sims)) sims, $(length(lifetimes)) inter-division times")

        if isempty(lifetimes)
            text!(ax_e, 0.5, 0.5; text = "no data",
                  align = (:center, :center), space = :relative,
                  fontsize = 9, color = :grey60)
            continue
        end

        hist!(ax_e, lifetimes; bins = 80, normalization = :pdf,
              color = (COL_DIV, 0.45), strokewidth = 0.3, strokecolor = COL_DIV)
        lines!(ax_e, collect(xs), pdf.(bd, xs);
               color = :black, linewidth = 1.5, linestyle = :dash)

        if d > 0.0
            text!(ax_e, 0.97, 0.95;
                text     = "death times not\nobservable (pruned)",
                align    = (:right, :top), space = :relative,
                fontsize = FS_ANNOT, color = COL_DEATH)
        end
    end

    linkyaxes!(axes_theory...);  linkxaxes!(axes_theory...)
    linkyaxes!(axes_empiric...); linkxaxes!(axes_empiric...)

    rowgap!(fig.layout, 1, 2.0)
    rowgap!(fig.layout, 2, 14.0)
    rowgap!(fig.layout, 3, 2.0)

    Legend(fig[5, 1:ncols],
        [
            LineElement(color = COL_DIV,   linewidth = 2.0),
            LineElement(color = COL_DEATH, linewidth = 2.0, linestyle = :dash),
            PolyElement(color = (COL_DIV, 0.45)),
            LineElement(color = :black,    linewidth = 1.5, linestyle = :dash),
        ],
        [label_div, label_death,
         "empirical inter-division times",
         "theoretical division PDF"];
        orientation  = :horizontal, labelsize    = FS_LEGEND,
        framevisible = false,       tellwidth    = false,
        patchsize    = (16.0, 8.0), nbanks       = 2)

    mkpath(dirname(outpath))
    save(outpath, fig)
    println("→ $outpath\n")
end

# ── Markov model ──────────────────────────────────────────────────────────────

println("\n── Markov timing distribution ─────────────────────────────────────")
make_timing_figure(
    dist_birth  = _ -> Exponential(1.0 / b),
    dist_death  = d -> Exponential(1.0 / d),
    raw_path    = d -> joinpath(RAW, "markov", "neutral_markov_N$(N_LOAD)_d$(d).jls"),
    label_div   = "division: Exponential(1/b)",
    label_death = "death: Exponential(1/d)  [theory only]",
    outpath     = joinpath(OUT, "markov", "timing_distributions.png"),
)

# ── Gamma model ───────────────────────────────────────────────────────────────

println("── Gamma timing distribution ──────────────────────────────────────")
make_timing_figure(
    dist_birth  = _ -> Gamma(k, 1.0 / (k * b)),
    dist_death  = d -> Gamma(k, 1.0 / (k * d)),
    raw_path    = d -> joinpath(RAW, "gamma", "neutral_gamma_N$(N_LOAD)_d$(d)_k$(k).jls"),
    label_div   = "division: Gamma($(Int(k)), 1/(k·b))",
    label_death = "death: Gamma($(Int(k)), 1/(k·d))  [theory only]",
    outpath     = joinpath(OUT, "gamma", "timing_distributions.png"),
)

println("All done.")
