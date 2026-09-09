using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using Distributions
using Test

include(joinpath(HELPERS, "divisiontimes.jl"))

# Unit suite for the observed inter-division time law derived in
# theory/divisiontimes.md.  Reads no simulation data: every target below is
# either a closed form from that note or a number published in its
# "Numerical verification" tables, so this runs in seconds and fails loudly if
# the quadrature or the Euler–Lotka bisection drifts.
#
#   g(τ)      = f_b(τ) (1 - F_d(τ))            defective "divides at age τ"
#   p_div     = ∫ g                            probability of dividing at all
#   realised  = g / p_div                      Effect 1 alone
#   observed  = 2 g(τ) exp(-λτ)                Effect 1 + 2, what a tree records
#   Euler–Lotka: 2 ∫ exp(-λτ) g(τ) dτ = 1      defines λ

const b = 1.0
const k = 5.0

@testset "observed inter-division times" begin

# ── Markov: everything is closed form ─────────────────────────────────────────
#
# g(τ) = b exp(-(b+d)τ), so p_div = b/(b+d), the realised lifetime is Exp(b+d),
# Euler–Lotka gives λ = b - d, and the tilt cancels d exactly: p_obs = Exp(2b)
# for every d.  theory/divisiontimes.md §"Markov (k = 1)".

@testset "markov d=$d" for d in (0.0, 0.5, 0.9)
    law = predict_observed_lifetimes(Exponential(1 / b),
                                     d == 0 ? nothing : Exponential(1 / d))

    @test law.λ     ≈ b - d          atol = 1e-8
    @test law.p_div ≈ b / (b + d)    atol = 1e-8

    # Realised lifetime is Exp(b+d) — the death clock speeds the observed clock.
    @test mean_realised(law) ≈ 1 / (b + d) atol = 1e-8
    for τ in (0.05, 0.5, 1.0, 3.0)
        @test pdf_realised(law, τ) ≈ pdf(Exponential(1 / (b + d)), τ) atol = 1e-10
    end

    # Observed lifetime is Exp(2b), independent of d — the factor-of-two the
    # figures show at τ = 0.
    @test mean_observed(law)   ≈ 1 / (2b) atol = 1e-8
    @test pdf_observed(law, 0.0) ≈ 2b     atol = 1e-8
    for τ in (0.05, 0.5, 1.0, 3.0)
        @test pdf_observed(law, τ) ≈ pdf(Exponential(1 / (2b)), τ) atol = 1e-10
    end
end

# ── Gamma k=5, d=0: closed form ───────────────────────────────────────────────
#
# λ = kb(2^{1/k} - 1) and p_obs = Gamma(k, 1/(k b 2^{1/k})): the tilt returns the
# same family with the same CV, rescaled by 2^{-1/k}.

@testset "gamma d=0 closed form" begin
    law = predict_observed_lifetimes(Gamma(k, 1 / (k * b)), nothing)

    @test law.λ     ≈ k * b * (2^(1 / k) - 1) atol = 1e-8
    @test law.λ     ≈ 0.74349                 atol = 1e-5   # published table
    @test law.p_div ≈ 1.0                     atol = 1e-10  # no death clock

    @test mean_observed(law) ≈ 2^(-1 / k) / b atol = 1e-8
    @test mean_observed(law) ≈ 0.87055        atol = 1e-5   # published table
    @test mean_realised(law) ≈ 1 / b          atol = 1e-8   # no death ⇒ untilted = f_b

    tilted = Gamma(k, 1 / (k * b * 2^(1 / k)))
    for τ in (0.1, 0.5, 0.87, 1.5, 3.0)
        @test pdf_observed(law, τ) ≈ pdf(tilted, τ)          atol = 1e-9
        @test pdf_realised(law, τ) ≈ pdf(Gamma(k, 1/(k*b)), τ) atol = 1e-9
    end
end

# ── Gamma k=5, d>0: no closed form, quadrature against the published table ────
#
# theory/divisiontimes.md §"Numerical verification", gamma rows.  1 - F_d is not
# exponential here, so g is not a Gamma and λ has to be found numerically.

@testset "gamma d=$d against published table" for (d, λ_pub, pdiv_pub, mean_pub, unt_pub) in
    ((0.5, 0.61226, 0.855, 0.83540, 0.92015),
     (0.9, 0.15614, 0.565, 0.76977, 0.78504))

    law = predict_observed_lifetimes(Gamma(k, 1 / (k * b)), Gamma(k, 1 / (k * d)))

    @test law.λ            ≈ λ_pub    atol = 1e-5
    @test law.p_div        ≈ pdiv_pub atol = 5e-4
    @test mean_observed(law) ≈ mean_pub atol = 1e-5
    @test mean_realised(law) ≈ unt_pub  atol = 1e-5
end

# ── Both densities are normalised ─────────────────────────────────────────────
#
# For p_obs this *is* Euler–Lotka: the equation that fixes λ is exactly the
# statement that 2 g exp(-λτ) integrates to 1.

@testset "normalisation" begin
    for (f_b, f_d) in ((Exponential(1 / b),   nothing),
                       (Exponential(1 / b),   Exponential(1 / 0.9)),
                       (Gamma(k, 1/(k*b)),    nothing),
                       (Gamma(k, 1/(k*b)),    Gamma(k, 1/(k*0.5))))
        law = predict_observed_lifetimes(f_b, f_d)
        @test integrate_density(τ -> pdf_observed(law, τ), law) ≈ 1.0 atol = 1e-8
        @test integrate_density(τ -> pdf_realised(law, τ), law) ≈ 1.0 atol = 1e-8
    end
end

# ── Deterministic: tilting a point mass leaves it where it is ─────────────────
#
# Euler–Lotka reads 2 exp(-λ/b) = 1, so λ = b log 2 — which is also the k → ∞
# limit of the gamma expression kb(2^{1/k} - 1).  All three laws coincide, which
# is why the deterministic panel shows no discrepancy at all.

@testset "deterministic" begin
    law = predict_observed_lifetimes(Dirac(1 / b), nothing)

    @test law.λ     ≈ b * log(2) atol = 1e-12
    @test law.p_div ≈ 1.0        atol = 1e-12

    @test mean_observed(law) ≈ 1 / b atol = 1e-12
    @test mean_realised(law) ≈ 1 / b atol = 1e-12

    # k → ∞ limit of the gamma case agrees.
    @test 500 * b * (2^(1 / 500) - 1) ≈ b * log(2) atol = 1e-3
end

# ── Guard rails ───────────────────────────────────────────────────────────────
#
# p_div ≤ 1/2 means fewer than one surviving daughter per cell on average: the
# process is not supercritical and Euler–Lotka has no positive root.  That is a
# modelling error, not a numerical one, so it must not return a silent zero.

@testset "non-supercritical input errors" begin
    # Markov with d > b: p_div = b/(b+d) = 1/3 < 1/2.
    @test_throws ErrorException predict_observed_lifetimes(
        Exponential(1 / b), Exponential(1 / (2b)))
end

end
