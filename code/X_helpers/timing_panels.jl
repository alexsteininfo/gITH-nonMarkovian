using CairoMakie
using Distributions
using Statistics

# Shared rendering for the "measured inter-division times" panels, used by both
# code/2_theory/neutral/plot_timing.jl (per-model, two-row figures) and
# plot_timing_combined.jl (all three models, measured row only).
#
# Kept out of plotting_functions.jl on purpose: that file is constants and
# computation helpers and deliberately does not depend on CairoMakie.
#
# Requires the caller to have included plotting_functions.jl (palette, font
# sizes) and divisiontimes.jl (the three laws).

# ── Curve colours for the three theoretical laws ──────────────────────────────
#
# Not added to the shared palette because plotting_functions.R mirrors it by hand
# and nothing on the R side draws these.

const COL_FB       = :black
const COL_REALISED = :grey45
const COL_POBS     = COL_THEORY_DET   # firebrick — the curve that should match

# ── Dirac (point-mass) rendering ──────────────────────────────────────────────
#
# The deterministic model draws its waiting times from Dirac(1/b), which has no
# density to plot: pdf() is 0 everywhere except a non-integrable spike at 1/b.
# Those panels therefore show *mass* rather than density, drawn as a stem with an
# arrowhead — the usual convention for δ(t − t₀).

function spike!(ax, t0, h; color, linewidth, linestyle = :solid, markersize = 16)
    lines!(ax, [t0, t0], [0.0, h]; color, linewidth, linestyle)
    scatter!(ax, [t0], [h]; color, marker = :utriangle, markersize)
end

commas(n::Integer) = replace(string(n), r"(?<=[0-9])(?=([0-9]{3})+$)" => ",")

# ── The measured panel ────────────────────────────────────────────────────────

"""
    draw_lifetime_panel!(ax, lifetimes, f_b, law, xs; annotate = true, bins = 80)

Histogram of measured inter-division times against all three theoretical laws
(theory/divisiontimes.md):

  f_b      the input division clock            black dashed
  g/p_div  Effect 1 — competing death clock    grey dash-dot
  p_obs    Effects 1+2 — plus the census tilt  firebrick solid

Only `p_obs` should track the histogram; the gap between the three is the point
of the panel.
"""
function draw_lifetime_panel!(ax, lifetimes, f_b, law, xs;
                              annotate::Bool = true, bins::Int = 80)
    hist!(ax, lifetimes; bins, normalization = :pdf,
          color = (COL_DIV, 0.45), strokewidth = 0.6, strokecolor = COL_DIV)
    lines!(ax, collect(xs), pdf.(f_b, xs);
           color = COL_FB, linewidth = 2.5, linestyle = :dash)
    lines!(ax, collect(xs), pdf_realised.(Ref(law), xs);
           color = COL_REALISED, linewidth = 2.5, linestyle = :dashdot)
    lines!(ax, collect(xs), pdf_observed.(Ref(law), xs);
           color = COL_POBS, linewidth = 3.0)

    # Measured against predicted, so the panel reads as a validation rather than
    # an assertion. f_b is listed last because it is the one the histogram is
    # *not* expected to match.
    annotate && text!(ax, 0.97, 0.95;
        text     = "λ = $(round(law.λ, digits = 4))\n" *
                   "mean measured   $(round(mean(lifetimes),    digits = 4))\n" *
                   "mean p_obs      $(round(mean_observed(law), digits = 4))\n" *
                   "mean g/p_div    $(round(mean_realised(law), digits = 4))\n" *
                   "mean f_b        $(round(mean(f_b),          digits = 4))",
        align    = (:right, :top), space = :relative,
        fontsize = FS_ANNOT, color = :grey40)
    return ax
end

"""
    draw_lifetime_spike_panel!(ax, lifetimes, f_b, law, xs; annotate = true)

The Dirac counterpart of [`draw_lifetime_panel!`](@ref). Every measured lifetime
should be exactly `t₀`, so a histogram would have no width; the measurement is
drawn as mass and the degeneracy is *checked* rather than assumed.

One stem, not three: the tilt exp(-λτ) cannot move a point mass, so f_b, g/p_div
and p_obs are the same Dirac here — which is why the deterministic model is the
one case with no discrepancy to explain.
"""
function draw_lifetime_spike_panel!(ax, lifetimes, f_b, law, xs; annotate::Bool = true)
    lo, hi = extrema(lifetimes)
    isapprox(lo, hi; atol = 1e-9) || error(
        "deterministic tree has non-constant inter-division times: " *
        "extrema = ($lo, $hi) — is this really the Dirac model?")
    t_obs = (lo + hi) / 2

    xlims!(ax, extrema(xs)...)
    ylims!(ax, 0.0, 1.3)
    # Wide translucent stem for the measurement, in place of the histogram bars
    # the other models get; theory dashed on top.
    spike!(ax, t_obs, 1.0; color = (COL_DIV, 0.45), linewidth = 10.0, markersize = 20)
    spike!(ax, f_b.value, 1.0;
           color = COL_FB, linewidth = 2.0, linestyle = :dash, markersize = 9)

    if annotate
        text!(ax, 0.97, 0.95;
            text     = "all $(commas(length(lifetimes))) inter-division\n" *
                       "times exactly t = $(t_obs)",
            align    = (:right, :top), space = :relative,
            fontsize = FS_ANNOT, color = :grey40)
        text!(ax, 0.97, 0.70;
            text     = "f_b = g/p_div = p_obs\n" *
                       "λ = $(round(law.λ, digits = 4)), but tilting a\n" *
                       "point mass leaves it in place",
            align    = (:right, :top), space = :relative,
            fontsize = FS_ANNOT, color = COL_POBS)
    end
    return ax
end

# ── Legend pieces ─────────────────────────────────────────────────────────────
#
# A vertical stem in the swatch for the point-mass panels, a horizontal line or
# patch for the density panels — matching what each actually draws.

stem_swatch(; kwargs...) =
    LineElement(; points = Point2f[(0.5, 0.0), (0.5, 1.0)], kwargs...)

# The three theoretical laws plus the histogram, in the order they are drawn.
lifetime_legend_entries() = (
    LegendElement[
        PolyElement(color = (COL_DIV, 0.45)),
        LineElement(color = COL_FB,       linewidth = 2.5, linestyle = :dash),
        LineElement(color = COL_REALISED, linewidth = 2.5, linestyle = :dashdot),
        LineElement(color = COL_POBS,     linewidth = 3.0),
    ],
    ["empirical inter-division times",
     "input clock  f_b",
     "realised  g/p_div  (death clock only)",
     "observed  p_obs = 2 g(t) exp(-λt)"],
)
