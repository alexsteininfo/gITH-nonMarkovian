using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using Distributions
using Serialization
using CairoMakie
using Statistics
using Random

include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "plotting_functions.jl"))
include(joinpath(HELPERS, "divisiontimes.jl"))
include(joinpath(HELPERS, "timing_panels.jl"))

# One figure holding only the *measured* row of each of the three timing models —
# the second row of figures/2_theory/neutral/{markov,gamma,deterministic}/
# timing_distributions.png, stacked. The per-model figures keep their theory row;
# this one drops it so the three models can be compared directly.
#
# Panels are drawn by the shared helpers in HELPERS/timing_panels.jl, so this
# figure and the per-model ones cannot drift apart.

# ── Parameters ────────────────────────────────────────────────────────────────

const b        = 1.0
const k        = 5.0
const D_VALUES = [0.0, 0.5, 0.9]
const N_LOAD     = 10_000
const N_LOAD_DET = 16_384   # deterministic sweep uses powers of two, not 10^n

const RAW = joinpath(DATA, "raw", "neutral")
const OUT = joinpath(FIGURES, "2_theory", "neutral", "combined", "timing", "timing.png")

xs = range(0.0, 6.0; length = 400)

# ── Rows ──────────────────────────────────────────────────────────────────────
#
# The deterministic model has no d > 0 case: its division time 1/b is always
# shorter than its death time 1/d, so no cell ever dies and the sweep was only
# ever run at d = 0. That row therefore has one panel and two empty cells rather
# than a fabricated comparison.

const ROWS = [
    (name       = "Markov  (k = 1)",
     dist_birth = _ -> Exponential(1.0 / b),
     dist_death = d -> Exponential(1.0 / d),
     raw_path   = d -> joinpath(RAW, "markov", "neutral_markov_N$(N_LOAD)_d$(d).jls"),
     d_values   = D_VALUES),

    (name       = "Gamma  (k = $(Int(k)))",
     dist_birth = _ -> Gamma(k, 1.0 / (k * b)),
     dist_death = d -> Gamma(k, 1.0 / (k * d)),
     raw_path   = d -> joinpath(RAW, "gamma",
                                "neutral_gamma_N$(N_LOAD)_d$(d)_k$(k).jls"),
     d_values   = D_VALUES),

    (name       = "Deterministic  (k = ∞)",
     dist_birth = _ -> Dirac(1.0 / b),
     dist_death = _ -> Dirac(Inf),          # never called: d = 0 only
     raw_path   = _ -> joinpath(RAW, "deterministic",
                                "neutral_deterministic_N$(N_LOAD_DET).jls"),
     d_values   = [0.0]),
]

# ── Figure ────────────────────────────────────────────────────────────────────

const NCOL = length(D_VALUES)
const NROW = length(ROWS)

# The deterministic row has only a d = 0 panel, so the bottom-most panel differs
# by column: column 1 ends in the deterministic row, columns 2–3 end in the gamma
# row. Label whichever panel is actually last, or columns 2–3 get no x label.
const LAST_ROW_IN_COL =
    [maximum(ri for ri in 1:NROW if ci <= length(ROWS[ri].d_values)) for ci in 1:NCOL]

fig = Figure(size = (NCOL * 430 + 150, NROW * 330 + 130), figure_padding = 18)

Label(fig[1, 2:(NCOL + 1)];
    text      = "Measured inter-division times vs. the three theoretical laws",
    fontsize  = FS_HEAD, font = :bold,
    tellwidth = false, padding = (0, 0, 4, 0))

# Density panels share a y-axis; the deterministic row shows probability *mass*
# on a 0–1.3 scale, so it is linked separately (i.e. not at all — one panel).
density_axes = Axis[]

for (ri, row) in enumerate(ROWS)
    dirac = row.dist_birth(first(row.d_values)) isa Dirac

    # Row label in its own narrow column, so it does not eat panel width.
    Label(fig[ri + 1, 1];
        text      = row.name,
        fontsize  = FS_TITLE, font = :bold, rotation = pi / 2,
        tellheight = false, padding = (0, 6, 0, 0))

    println("── $(row.name) ─────────────────────────────────")

    for (ci, d) in enumerate(row.d_values)
        bd  = row.dist_birth(d)
        law = predict_observed_lifetimes(bd, d > 0.0 ? row.dist_death(d) : nothing)

        ax = Axis(fig[ri + 1, ci + 1];
            title          = ri == 1 ? "d = $d" : "",
            titlesize      = FS_TITLE,
            xlabel         = ri == LAST_ROW_IN_COL[ci] ? "inter-division time" : "",
            ylabel         = ci > 1  ? ""                :
                             dirac   ? "probability mass" : "probability density",
            xlabelsize     = FS_LABEL, ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,  yticklabelsize = FS_TICK,
            spinewidth     = 1.3)

        # The deterministic row sits under column 1, which is titled "d = 0.0"
        # only on the top row; restate it so the panel is readable on its own.
        ri > 1 && (ax.title = "d = $d")
        ri > 1 && (ax.titlesize = FS_TICK + 2)

        fpath = row.raw_path(d)
        if !isfile(fpath)
            text!(ax, 0.5, 0.5; text = "data not found",
                  align = (:center, :center), space = :relative,
                  fontsize = FS_LABEL, color = :grey60)
            continue
        end

        sims      = deserialize(fpath)
        lifetimes = vcat([celllifetimes(something(s.tree_root))
                          for s in sims if !isnothing(s.tree_root)]...)
        println("  d=$d: $(length(sims)) sims, $(length(lifetimes)) inter-division times")

        if isempty(lifetimes)
            text!(ax, 0.5, 0.5; text = "no data",
                  align = (:center, :center), space = :relative,
                  fontsize = FS_LABEL, color = :grey60)
            continue
        end

        if dirac
            draw_lifetime_spike_panel!(ax, lifetimes, bd, law, xs)
        else
            draw_lifetime_panel!(ax, lifetimes, bd, law, xs)
            push!(density_axes, ax)
        end
    end

    # Deterministic: say why the other two columns are absent rather than
    # leaving the reader to wonder whether the run is missing.
    if length(row.d_values) < NCOL
        Label(fig[ri + 1, (length(row.d_values) + 2):(NCOL + 1)];
            text      = "no d > 0 panels: with Dirac timing the division time 1/b\n" *
                        "always precedes the death time 1/d, so no cell ever dies",
            fontsize  = FS_ANNOT, color = :grey50,
            tellwidth = false, tellheight = false)
    end
    println()
end

linkxaxes!(density_axes...)
linkyaxes!(density_axes...)

leg_elems, leg_labels = lifetime_legend_entries()
Legend(fig[NROW + 2, 1:(NCOL + 1)], leg_elems, leg_labels;
    orientation  = :horizontal, labelsize    = FS_LEGEND,
    framevisible = false,       tellwidth    = false,
    patchsize    = (28.0, 14.0), nbanks      = 1)

colgap!(fig.layout, 1, 2.0)
rowgap!(fig.layout, 1, 4.0)
rowgap!(fig.layout, NROW + 1, 6.0)

# The deterministic row carries one spike and two empty cells; at full height it
# reads as a layout error rather than a real asymmetry of the model.
rowsize!(fig.layout, NROW + 1, Auto(0.62))

mkpath(dirname(OUT))
save(OUT, fig)
println("→ $OUT")
