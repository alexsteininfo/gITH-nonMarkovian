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

# ── Parameters ────────────────────────────────────────────────────────────────

const b        = 1.0
const k        = 5.0
const D_VALUES = [0.0, 0.5, 0.9]
const N_LOAD   = 10_000   # use largest trees for richest empirical distribution
const N_LOAD_DET = 16_384 # deterministic sweep uses powers of two, not 10^n

const RAW = joinpath(DATA,    "raw", "neutral")
const OUT = joinpath(FIGURES, "2_theory", "neutral")

xs = range(0.0, 6.0; length = 400)

# Panel rendering (spike!, commas, the three-law overlay, the legend swatches)
# lives in HELPERS/timing_panels.jl, shared with plot_timing_combined.jl.

# ── Figure factory ────────────────────────────────────────────────────────────
#
# Produces a two-row timing figure:
#   Top row    — theoretical PDFs for birth (and death for d > 0), i.e. the
#                clocks the simulation was *given*
#   Bottom row — empirical inter-division times from simulation trees, against
#                all three theoretical laws (see theory/divisiontimes.md):
#                  f_b       the input division clock            black dashed
#                  g/p_div   Effect 1 — competing death clock    grey dash-dot
#                  p_obs     Effects 1+2 — plus the census tilt  firebrick solid
#                Only p_obs should track the histogram; the gap between the
#                three is the point of the row.
#
# dist_birth(d)  : returns the birth distribution for a given d
# dist_death(d)  : returns the death distribution (only called for d > 0)
# raw_path(d)    : returns the raw .jls file path for a given d
# label_div      : legend label for the division distribution line
# label_death    : legend label for the death distribution line
# outpath        : where to save the figure
# d_values       : death rates to show, one column each
#
# A Dirac birth distribution (the deterministic model) switches both rows to the
# point-mass rendering described above; d_values is then [0.0] and no death
# distribution is drawn, so the death legend entry is dropped.

function make_timing_figure(;
    dist_birth,
    dist_death,
    raw_path,
    label_div,
    label_death,
    outpath,
    d_values = D_VALUES,
)
    ncols = length(d_values)
    dirac = dist_birth(first(d_values)) isa Dirac
    fig   = Figure(size = (max(ncols, 2) * 380, 2 * 320 + 110), figure_padding = 18)

    Label(fig[1, 1:ncols];
        text      = dirac ? "Implemented (theoretical timing distribution)" :
                            "Implemented (theoretical PDFs)",
        fontsize  = FS_HEAD, font = :bold,
        tellwidth = false, padding = (0, 0, 4, 0))
    Label(fig[3, 1:ncols];
        text      = "Measured (inter-division times from simulation trees)",
        fontsize  = FS_HEAD, font = :bold,
        tellwidth = false, padding = (0, 8, 4, 0))

    axes_theory  = Vector{Axis}(undef, ncols)
    axes_empiric = Vector{Axis}(undef, ncols)

    for (ci, d) in enumerate(d_values)
        bd = dist_birth(d)

        # The two derived laws. Solving Euler–Lotka is a quadrature plus a
        # bisection, so it costs well under a second per column — negligible
        # beside deserializing the trees below.
        law = predict_observed_lifetimes(bd, d > 0.0 ? dist_death(d) : nothing)

        # Row 1: theoretical ─────────────────────────────────────────────────
        ax_t = Axis(fig[2, ci];
            title          = "d = $d",
            titlesize      = FS_TITLE,
            xlabel         = "waiting time",
            ylabel         = ci > 1        ? ""                 :
                             dirac         ? "probability mass"  : "probability density",
            xlabelsize     = FS_LABEL, ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,  yticklabelsize = FS_TICK,
            spinewidth     = 1.3)
        axes_theory[ci] = ax_t

        if dirac
            # A lone vertical line has no x-extent to autoscale from, so the
            # window is pinned to the same range the other models are drawn on.
            xlims!(ax_t, extrema(xs)...)
            ylims!(ax_t, 0.0, 1.3)
            spike!(ax_t, bd.value, 1.0; color = COL_DIV, linewidth = 4.0)
        else
            lines!(ax_t, collect(xs), pdf.(bd, xs); color = COL_DIV, linewidth = 3.0)
        end

        if d > 0.0
            dd = dist_death(d)
            lines!(ax_t, collect(xs), pdf.(dd, xs);
                   color = COL_DEATH, linewidth = 3.0, linestyle = :dash)
        end

        # Row 2: empirical ───────────────────────────────────────────────────
        ax_e = Axis(fig[4, ci];
            title          = "d = $d",
            titlesize      = FS_TITLE,
            xlabel         = "inter-division time",
            ylabel         = ci > 1        ? ""                 :
                             dirac         ? "probability mass"  : "probability density",
            xlabelsize     = FS_LABEL, ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,  yticklabelsize = FS_TICK,
            spinewidth     = 1.3)
        axes_empiric[ci] = ax_e

        fpath = raw_path(d)
        if !isfile(fpath)
            text!(ax_e, 0.5, 0.5; text = "data not found",
                  align = (:center, :center), space = :relative,
                  fontsize = FS_LABEL, color = :grey60)
            continue
        end

        sims      = deserialize(fpath)
        lifetimes = vcat([celllifetimes(something(s.tree_root))
                          for s in sims if !isnothing(s.tree_root)]...)
        println("  d=$d: $(length(sims)) sims, $(length(lifetimes)) inter-division times")

        if isempty(lifetimes)
            text!(ax_e, 0.5, 0.5; text = "no data",
                  align = (:center, :center), space = :relative,
                  fontsize = FS_LABEL, color = :grey60)
            continue
        end

        if dirac
            draw_lifetime_spike_panel!(ax_e, lifetimes, bd, law, xs)
        else
            draw_lifetime_panel!(ax_e, lifetimes, bd, law, xs)
        end
    end

    linkyaxes!(axes_theory...);  linkxaxes!(axes_theory...)
    linkyaxes!(axes_empiric...); linkxaxes!(axes_empiric...)

    rowgap!(fig.layout, 1, 2.0)
    rowgap!(fig.layout, 2, 14.0)
    rowgap!(fig.layout, 3, 2.0)

    leg_elems, leg_labels = if dirac
        (LegendElement[
             stem_swatch(color = COL_DIV, linewidth = 4.0),
             stem_swatch(color = (COL_DIV, 0.45), linewidth = 10.0),
             stem_swatch(color = COL_FB,  linewidth = 2.0, linestyle = :dash),
         ],
         [label_div,
          "empirical inter-division times",
          "f_b = g/p_div = p_obs  (all three coincide)"])
    else
        lt_elems, lt_labels = lifetime_legend_entries()
        (LegendElement[
             LineElement(color = COL_DIV,   linewidth = 3.0),
             LineElement(color = COL_DEATH, linewidth = 3.0, linestyle = :dash),
             lt_elems...,
         ],
         [label_div, label_death, lt_labels...])
    end

    Legend(fig[5, 1:ncols], leg_elems, leg_labels;
        orientation  = :horizontal, labelsize    = FS_LEGEND,
        framevisible = false,       tellwidth    = false,
        patchsize    = (28.0, 14.0), nbanks       = 2)

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

# ── Deterministic model ───────────────────────────────────────────────────────
#
# Dirac(1/b) division timing and no death — the same distributions the
# simulation uses (see 1_data_generation/1_simulation_runs/neutral/
# growth_neutral_deterministic.jl), so d = 0 is the only column.

println("── Deterministic timing distribution ──────────────────────────────")
make_timing_figure(
    dist_birth  = _ -> Dirac(1.0 / b),
    dist_death  = _ -> Dirac(Inf),           # never called: d = 0 only
    raw_path    = _ -> joinpath(RAW, "deterministic",
                                "neutral_deterministic_N$(N_LOAD_DET).jls"),
    label_div   = "division: Dirac(1/b)",
    label_death = "",                        # no death entry in the legend
    outpath     = joinpath(OUT, "deterministic", "timing_distributions.png"),
    d_values    = [0.0],
)

println("All done.")
