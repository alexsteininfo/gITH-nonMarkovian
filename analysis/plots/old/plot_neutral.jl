using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Distributions
using Serialization
using CairoMakie
using Statistics

include(joinpath(dirname(@__DIR__), "helpers", "types.jl"))

# ── Parameters ────────────────────────────────────────────────────────────────

const b        = 1.0
const k        = 5.0
const D_VALUES = [0.0, 0.5, 0.9]
const N_LOAD   = 10_000   # use the largest trees for empirical distributions

const RAW_DIR = joinpath(@__DIR__, "..", "..", "data", "raw", "growth_neutral")
const OUTDIR  = joinpath(@__DIR__, "..", "..", "plots", "neutral")

const COL_DIV    = :royalblue
const COL_DEATH  = :darkorange
const COL_THEORY = :black

const FS_HEAD   = 11
const FS_TITLE  = 10
const FS_LABEL  =  9
const FS_TICK   =  8
const FS_LEGEND =  8
const FS_ANNOT  =  7

mkpath(OUTDIR)

ncols = length(D_VALUES)   # 3
nrows = 2                  # implemented / measured

fig = Figure(size = (ncols * 300, nrows * 240 + 80), figure_padding = 14)

Label(fig[1, 1:ncols];
    text      = "Implemented (theoretical Gamma PDFs)",
    fontsize  = FS_HEAD,
    font      = :bold,
    tellwidth = false,
    padding   = (0, 0, 4, 0),
)
Label(fig[3, 1:ncols];
    text      = "Measured (inter-division times from simulation trees)",
    fontsize  = FS_HEAD,
    font      = :bold,
    tellwidth = false,
    padding   = (0, 8, 4, 0),
)

xs = range(0.0, 6.0; length = 400)
div_dist = Gamma(k, 1.0 / (k * b))

axes_theory  = Vector{Axis}(undef, ncols)
axes_empiric = Vector{Axis}(undef, ncols)

for (ci, d) in enumerate(D_VALUES)
    # ── Row 1: theoretical ───────────────────────────────────────────────────

    ax_t = Axis(fig[2, ci];
        title          = "d = $d",
        titlesize      = FS_TITLE,
        xlabel         = "waiting time",
        ylabel         = ci == 1 ? "probability density" : "",
        xlabelsize     = FS_LABEL,
        ylabelsize     = FS_LABEL,
        xticklabelsize = FS_TICK,
        yticklabelsize = FS_TICK,
        spinewidth     = 0.8,
    )
    axes_theory[ci] = ax_t

    lines!(ax_t, collect(xs), pdf.(div_dist, xs);
           color = COL_DIV, linewidth = 2.0, label = "division  Gamma($k, $(round(1/(k*b), digits=2)))")

    if d > 0.0
        death_dist = Gamma(k, 1.0 / (k * d))
        lines!(ax_t, collect(xs), pdf.(death_dist, xs);
               color = COL_DEATH, linewidth = 2.0, linestyle = :dash,
               label = "death  Gamma($k, $(round(1/(k*d), digits=2)))")
    end

    # ── Row 2: empirical ─────────────────────────────────────────────────────

    ax_e = Axis(fig[4, ci];
        title          = "d = $d",
        titlesize      = FS_TITLE,
        xlabel         = "inter-division time",
        ylabel         = ci == 1 ? "probability density" : "",
        xlabelsize     = FS_LABEL,
        ylabelsize     = FS_LABEL,
        xticklabelsize = FS_TICK,
        yticklabelsize = FS_TICK,
        spinewidth     = 0.8,
    )
    axes_empiric[ci] = ax_e

    fpath = joinpath(RAW_DIR, "growth_neutral_N$(N_LOAD)_d$(d)_k$(k).jls")
    if !isfile(fpath)
        text!(ax_e, 0.5, 0.5; text = "data not found", align = (:center, :center),
              space = :relative, fontsize = 9, color = :grey60)
        continue
    end

    sims = deserialize(fpath)
    println("d=$d: $(length(sims)) sims loaded")

    lifetimes = vcat([celllifetimes(something(s.tree_root))
                      for s in sims if !isnothing(s.tree_root)]...)
    println("  → $(length(lifetimes)) inter-division times")

    if isempty(lifetimes)
        text!(ax_e, 0.5, 0.5; text = "no data", align = (:center, :center),
              space = :relative, fontsize = 9, color = :grey60)
        continue
    end

    hist!(ax_e, lifetimes; bins = 80, normalization = :pdf,
          color = (COL_DIV, 0.45), strokewidth = 0.3, strokecolor = COL_DIV)
    lines!(ax_e, collect(xs), pdf.(div_dist, xs);
           color = COL_THEORY, linewidth = 1.5, linestyle = :dash)

    if d > 0.0
        text!(ax_e, 0.97, 0.95;
            text      = "death times not\nobservable (pruned)",
            align     = (:right, :top),
            space     = :relative,
            fontsize  = FS_ANNOT,
            color     = COL_DEATH,
        )
    end
end

# ── Link axes ─────────────────────────────────────────────────────────────────

linkyaxes!(axes_theory...)
linkxaxes!(axes_theory...)
linkyaxes!(axes_empiric...)
linkxaxes!(axes_empiric...)

# ── Row gaps ─────────────────────────────────────────────────────────────────

rowgap!(fig.layout, 1, 2.0)
rowgap!(fig.layout, 2, 14.0)
rowgap!(fig.layout, 3, 2.0)

# ── Legend ───────────────────────────────────────────────────────────────────

Legend(fig[5, 1:ncols],
    [
        LineElement(color = COL_DIV,   linewidth = 2.0),
        LineElement(color = COL_DEATH, linewidth = 2.0, linestyle = :dash),
        PolyElement(color = (COL_DIV, 0.45)),
        LineElement(color = COL_THEORY, linewidth = 1.5, linestyle = :dash),
    ],
    [
        "division: Gamma($k, 1/(k·b))",
        "death: Gamma($k, 1/(k·d))  [theory only]",
        "empirical inter-division times",
        "theoretical division PDF",
    ];
    orientation  = :horizontal,
    labelsize    = FS_LEGEND,
    framevisible = false,
    tellwidth    = false,
    patchsize    = (16.0, 8.0),
    nbanks       = 2,
)

# ── Save ─────────────────────────────────────────────────────────────────────

outfile = joinpath(OUTDIR, "timing_distributions.png")
save(outfile, fig)
println("→ $outfile")
