using Distributions

# Observed inter-division times in a growing tree — theory/divisiontimes.md.
#
# The lifetimes a lineage tree records are NOT draws from the division clock the
# simulation was given.  Two filters sit between them:
#
#   Effect 1  a cell only divides if its division clock beats its death clock,
#             giving the defective density  g(τ) = f_b(τ) (1 - F_d(τ)),
#             normalised to  g/p_div  with  p_div = ∫ g ≤ 1;
#   Effect 2  a cell only *contributes* a lifetime if it finished dividing before
#             the census.  Births accumulate like exp(λs) in a growing
#             population, so the eligible birth window tilts the density by
#             exp(-λτ) — the slow cyclers are the ones caught unfinished.
#
# Together:   p_obs(τ) = 2 g(τ) exp(-λτ),   with λ from Euler–Lotka
#                        2 ∫ exp(-λτ) g(τ) dτ = 1.
#
# The Euler–Lotka condition both defines λ and is exactly the statement that
# p_obs normalises, which is where the leading 2 comes from.
#
# Consequence worth repeating: do not fit a clock to celllifetimes() output.
# Fitting Exponential to Markov trees returns 2b; fitting Gamma(k,·) to gamma
# trees returns the right k but 2^{1/k} b.
#
# There is deliberately ONE numerical path here.  The Markov and gamma-d=0 cases
# have closed forms (λ = b - d with p_obs = Exp(2b); λ = kb(2^{1/k} - 1) with
# p_obs = Gamma(k, 1/(k b 2^{1/k}))), but short-circuiting them in production
# would leave the general quadrature checked only against the 5-decimal table in
# the note.  They are used as test oracles in verify_divisiontimes.jl instead, so
# the code path the gamma-d>0 panels rely on is the one the exact cases verify.

# ── Quadrature ────────────────────────────────────────────────────────────────
#
# Composite Simpson on [0, τ_max].  The integrands are smooth and decay at least
# as fast as f_b, so a fixed grid is both adequate and reproducible; τ_max comes
# from an extreme quantile of the division clock rather than a hard-coded bound,
# so the same defaults work for Exponential(1) and Gamma(5, 1/5) alike.

const DIVTIME_PANELS = 20_000   # Simpson panels (even); h ≈ 2.6e-3 for Exp(1)
const DIVTIME_TAILQ  = 1e-15    # τ_max = 1.5 · quantile(f_b, 1 - DIVTIME_TAILQ)

function simpson(f, a::Float64, b::Float64, n::Int)
    iseven(n) || (n += 1)
    h = (b - a) / n
    s = f(a) + f(b)
    for i in 1:(n - 1)
        s += (isodd(i) ? 4.0 : 2.0) * f(a + i * h)
    end
    return s * h / 3
end

# ── The law ───────────────────────────────────────────────────────────────────
#
# f_d === nothing means no death clock (d = 0), which is distinct from a death
# clock that never fires: it keeps 1 - F_d ≡ 1 exact rather than relying on
# Exponential(Inf) behaving.

struct ObservedLifetimeLaw{B,D}
    f_b::B
    f_d::D
    λ::Float64
    p_div::Float64
    τ_max::Float64
    panels::Int
end

isdirac(law::ObservedLifetimeLaw) = law.f_b isa Dirac

# Defective density of "divides at age τ", g(τ) = f_b(τ)(1 - F_d(τ)).
function defective_density(f_b, f_d, τ::Float64)
    p = pdf(f_b, τ)
    p == 0.0 && return 0.0                     # short-circuits the Gamma tail
    return f_d === nothing ? p : p * ccdf(f_d, τ)
end

"""
    predict_observed_lifetimes(f_b, f_d = nothing)

Solve Euler–Lotka for the population growth rate `λ` and return the law of the
inter-division times a lineage tree records, given division clock `f_b` and
death clock `f_d` (`nothing` for no death).

Query it with [`pdf_realised`](@ref) (Effect 1 alone: the lifetime of a cell that
does divide) and [`pdf_observed`](@ref) (Effect 1 + 2: what the tree records),
and their `mean_realised` / `mean_observed` counterparts.
"""
function predict_observed_lifetimes(f_b, f_d = nothing;
                                    panels::Int = DIVTIME_PANELS,
                                    tailq::Float64 = DIVTIME_TAILQ)

    # Dirac division clock: g is a point mass, so there is nothing to integrate.
    # Euler–Lotka collapses to 2 p_div exp(-λ t₀) = 1.  Tilting a point mass
    # leaves it in place, so all three laws coincide — the reason the
    # deterministic panel shows no discrepancy.
    if f_b isa Dirac
        t0 = float(f_b.value)
        t0 > 0 || error("Dirac division clock at t₀ = $t0; expected t₀ > 0.")
        p_div = f_d === nothing ? 1.0 : float(ccdf(f_d, t0))
        check_supercritical(p_div)
        return ObservedLifetimeLaw(f_b, f_d, log(2 * p_div) / t0, p_div, t0, 0)
    end

    τ_max = 1.5 * float(quantile(f_b, 1 - tailq))
    g(τ)  = defective_density(f_b, f_d, τ)

    p_div = simpson(g, 0.0, τ_max, panels)
    check_supercritical(p_div)

    # Φ(λ) = 2 ∫ exp(-λτ) g(τ) dτ - 1 is strictly decreasing, Φ(0) = 2p_div - 1 > 0.
    Φ(λ) = 2 * simpson(τ -> exp(-λ * τ) * g(τ), 0.0, τ_max, panels) - 1

    lo, hi = 0.0, 1.0
    while Φ(hi) > 0
        hi *= 2
        hi > 1e12 && error("Euler–Lotka bisection failed to bracket λ below 1e12.")
    end
    for _ in 1:200                       # ~1e-13 relative on any plausible λ
        mid = (lo + hi) / 2
        Φ(mid) > 0 ? (lo = mid) : (hi = mid)
    end

    return ObservedLifetimeLaw(f_b, f_d, (lo + hi) / 2, p_div, τ_max, panels)
end

# p_div ≤ 1/2 means fewer than one surviving daughter per division on average:
# the process is not supercritical, no positive λ solves Euler–Lotka, and the
# "growing population" the tilt describes does not exist.  A modelling error
# rather than a numerical one, so it must not return a silent zero.
function check_supercritical(p_div::Float64)
    p_div > 0.5 || error(
        "p_div = $(round(p_div, digits = 6)) ≤ 1/2: the process is not " *
        "supercritical, so Euler–Lotka has no positive root. Check that the " *
        "death clock is slower than the division clock.")
end

dirac_guard(law) = isdirac(law) && error(
    "$(law.f_b) has no density — its observed law is a point mass at " *
    "$(law.f_b.value). Plot it as mass, or use mean_observed/mean_realised.")

"""
    pdf_realised(law, τ)

Effect 1 alone: `g(τ)/p_div`, the age at which a cell divides *given that it
divides at all*.  Still a per-cell property — this is what a single tracked
lineage would report, with no census bias.
"""
function pdf_realised(law::ObservedLifetimeLaw, τ::Real)
    dirac_guard(law)
    return defective_density(law.f_b, law.f_d, float(τ)) / law.p_div
end

"""
    pdf_observed(law, τ)

Effects 1 and 2: `p_obs(τ) = 2 g(τ) exp(-λτ)`, the density of a division event
drawn uniformly from those a tree censused at fixed `N` actually records.  This
is the curve a measured `celllifetimes` histogram converges to.
"""
function pdf_observed(law::ObservedLifetimeLaw, τ::Real)
    dirac_guard(law)
    t = float(τ)
    return 2 * defective_density(law.f_b, law.f_d, t) * exp(-law.λ * t)
end

"""
    integrate_density(f, law)

Integrate `f` over the law's own quadrature grid — used to confirm that both
densities normalise (for `pdf_observed` that *is* Euler–Lotka) and to take means.
"""
function integrate_density(f, law::ObservedLifetimeLaw)
    dirac_guard(law)
    return simpson(τ -> f(τ), 0.0, law.τ_max, law.panels)
end

"""
    mean_realised(law)

Mean of `g/p_div`.  For Markov this is `1/(b+d)`; the death clock speeds the
realised clock without touching the tilt.
"""
mean_realised(law::ObservedLifetimeLaw) =
    isdirac(law) ? float(law.f_b.value) :
    integrate_density(τ -> τ * pdf_realised(law, τ), law)

"""
    mean_observed(law)

Mean of `p_obs`.  For Markov this is `1/(2b)` for every `d`; for gamma with
`d = 0` it is `2^{-1/k}/b`, the factor that explains the whole figure family.
"""
mean_observed(law::ObservedLifetimeLaw) =
    isdirac(law) ? float(law.f_b.value) :
    integrate_density(τ -> τ * pdf_observed(law, τ), law)
