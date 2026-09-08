using Distributions

# Theorem 3 from: Genetics 230(4) iyaf101 (2025)
# E[M_j] = N * sum_l Poisson(l; averageGen) * Poisson(j; l*m)

# ── SFS and TMB (Eq. 9 / Gunnarsson, Stein & Werner 2025) ────────────────────
#
# omega = 2*mu throughout, where mu is the per-division mutation rate.
# The factor of 2 arises because the formula counts mutations on ALL branches
# of the phylogenetic tree — both internal branches (shared mutations, k≥2) and
# leaf branches (private mutations, k=1).  sitefrequencyspectrum fills sfs[1] from the
# leaf branches too, so empirical and theoretical spectra are directly
# comparable at every k, including k=1.
#
# Derivation: for d=0, E[sfs[k]] ≈ 2·nu·N/(k(k+1)) empirically (factor-of-2
# over the formula's N/(k(k+1))), confirmed by pooling 200 simulations.

# SFS for a perfect binary tree (deterministic Dirac timing, N must be a power of 2).
# At depth d from the root there are 2^d branches, each subtending N/2^d leaves.
# With ν mutations per division: E[SFS[k]] = (N/k)·ν for k ∈ {1,2,4,…,N/2}; 0 otherwise.
function predict_sfs_deterministic(k::Int, N::Int, nu::Float64)
    if k < 1 || k > N ÷ 2 || (k & (k - 1)) != 0
        return 0.0
    end
    return N * nu / k
end

# Expected SFS at frequency k (cells), Markovian birth-death process.
# Returns expected number of distinct mutations found in exactly k of the N
# final cells (Eq. 9, Stein & Werner 2025).
function predict_sfs_theory(k::Int, b::Float64, d::Float64, N::Int, mu::Float64)
    omega = 2.0 * mu
    p0 = d / b
    if p0 == 0.0
        return omega * N / (k * (k + 1))
    end
    total = 0.0
    for j in 0:1000
        add = p0^j / ((k + j) * (k + j + 1))
        total += add
        add <= 1e-10 && break
    end
    return omega * N * total
end

# Expected total number of distinct mutations in the final population,
# Markovian birth-death process (Gunnarsson's result; Stein & Werner 2025).
# Counts mutations on all branches (internal + leaf): TMB = sum_k E[S_k].
function predict_tmb_theory(b::Float64, d::Float64, N::Int, mu::Float64)
    omega = 2.0 * mu
    if d == 0.0
        return omega * (N - 1)
    else
        rho   = d / b
        sigma = 1.0 - rho
        MB    = -log(sigma - rho / N) / rho
        return omega * N * MB
    end
end

function predict_averageGeneration(b, d, N)
    γ = 0.5772156649015329   # Euler–Mascheroni constant
    if d == 0
        return 2 * (γ + log(N) - 1)
    else
        tN = (γ + log(N * (b - d) / b)) / (b - d)
        return 2 * b * tN - 2
    end
end

function predict_scMBdist(j_axis, averageGen, m, N)
    Lmax  = round(Int, 4 * averageGen)
    D_div = Poisson(averageGen)
    return Float64[
        N * sum(pdf(D_div, l) * pdf(Poisson(m * l), j) for l in 1:Lmax)
        for j in j_axis
    ]
end

# Deterministic mean divisional depth for the Dirac timing model: the population
# doubles each generation (synchronous binary tree), so depth = log₂(N) = log(N)/log(2).
# This holds for any d < b: since division time 1/b < death time 1/d, division always
# precedes death deterministically and no deaths occur.  The actual population overshoots
# to the next power of 2 above N_target; pass actual_N (= 2^k) for an exact result.
function predict_averageGeneration_deterministic(N)
    log(N) / log(2)
end

# Deterministic scMB: divisional distribution is a point mass → equation 24
# collapses to a single Poisson(m·⟨l⟩_det).
function predict_scMBdist_deterministic(j_axis, avgen_det, m, N)
    D = Poisson(m * avgen_det)
    return Float64[N * pdf(D, j) for j in j_axis]
end
