# Theoretical results used in this project

Reference: Stein & Werner (2025), *Genetics* **230**(4): iyaf101.
https://doi.org/10.1093/genetics/iyaf101

---

## Notation

| Symbol | Meaning |
|--------|---------|
| N | Final population size |
| b | Birth (division) rate |
| d | Death rate |
| m | Mean neutral mutations per division (Poisson) |
| l | Number of divisions a cell has undergone (divisional depth) |
| D_l | Number of cells with exactly l divisions |
| M_j | Number of cells carrying exactly j mutations |
| ⟨l⟩ | Mean divisional depth across all cells |

---

## Theorem 3 — Single-cell mutational burden distribution (general)

For a population of N cells, the expected scMB distribution is expressed in terms
of the expected divisional distribution {E[D_l]} by

$$E[M_j] = \sum_{l=1}^{\infty} E[D_l] \cdot \frac{(lm)^j \, e^{-lm}}{j!}
         = \sum_{l=1}^{\infty} E[D_l] \cdot \text{Pois}(j;\, lm)$$

**(Eq. 24, Stein & Werner 2025)**

**Interpretation.** Each cell that has undergone l divisions accumulates mutations
independently across divisions.  Because each division contributes Poisson(m)
neutral mutations, l divisions yield Poisson(l·m) total mutations.  The
population-level distribution is therefore a mixture of Poisson distributions
weighted by the divisional distribution.

Theorem 3 is **model-agnostic**: the formula holds for any birth–death process;
only the divisional distribution E[D_l] depends on the timing model.

---

## Markovian case (exponential / Gillespie)

When division and death times are exponentially distributed, the divisional
distribution at the time the population first reaches size N follows a Poisson
distribution with corrected mean ⟨l⟩ **(mean-field approximation, Eq. 25,
Stein & Werner 2025)**:

$$P_l = E[D_l] / N \sim \text{Pois}(\langle l \rangle), \qquad
  \langle l \rangle = 2b\,t_N - 2$$

where t_N is the mean time to reach size N (Eq. 26, Stein & Werner 2025).

Plugging into Theorem 3 gives a **compound Poisson** distribution:

$$E[M_j] = N \sum_{l=1}^{\infty} \text{Pois}(l;\,\langle l \rangle)
            \cdot \text{Pois}(j;\, lm)$$

This mixture is overdispersed relative to a plain Poisson.  By the law of total
variance, with L ~ Pois(⟨l⟩):

$$\text{Var}[M] = m\langle l\rangle + m^2\langle l\rangle
               = m\langle l\rangle(1 + m)$$

$$\text{Index of dispersion} = \frac{\text{Var}[M]}{E[M]} = 1 + m$$

The mutation-rate estimator derived from this is

$$\hat{m} = \frac{\text{Var}[\text{mpc}]}{E[\text{mpc}]} - 1$$

which is unbiased for exponential trees and systematically underestimates m for
more balanced (non-Markovian) trees, since σ²_l < ⟨l⟩ there.

### Limitation: mean-field approximation breaks down for d > 0

Eq. 25 (L ~ Poisson) is a mean-field approximation proven in the appendix of
Stein & Werner (2025).  The paper itself notes that the prediction is imperfect
when deaths are present.  The approximation assumes each cell independently
samples its divisional depth from the population mean, **ignoring correlations
between cells that share common ancestors**.

**Mechanism.**  Every division event on the path from the root to any cell is
shared by that cell and all its descendants.  Early division events — occurring
when the population is still small — are therefore shared by the largest numbers
of descendants and create the strongest positive correlations among depths.  In
the Yule process (d = 0) this effect is mild because all lineages grow
symmetrically and no single branch dominates.  When d > 0, deaths prune many
branches, and the surviving cells tend to coalesce onto a few early ancestors
whose division history is now inherited by a large fraction of the final
population.  This creates highly heterogeneous depth distributions within a
single tree: "ancient" lineages are deep; recently established lineages are
shallow.  The mean-field approximation, which treats each cell's depth as an
independent Poisson draw, cannot capture this heterogeneity.

**Empirical evidence.**  The index of dispersion D_L = Var[L] / E[L] of
empirical leaf depths (pooled over 200 simulations each):

| N | d | ⟨L⟩ predicted | ⟨L⟩ empirical | D_L = Var[L]/E[L] |
|---|---|---|---|---|
| 1 000 | 0.0 | 12.97 | 12.95 | 0.79 |
| 1 000 | 0.5 | 25.17 | 25.04 | 1.22 |
| 1 000 | 0.9 | 101.65 | 101.91 | **7.23** |
| 10 000 | 0.0 | 17.58 | 17.56 | 0.85 |
| 10 000 | 0.5 | 34.38 | 34.27 | 1.16 |
| 10 000 | 0.9 | 147.70 | 148.08 | **5.39** |

The mean prediction is accurate for all d; the variance is
underestimated at high d.  For d = 0.9 and μ = 2, the actual scMB index of
dispersion D_M ≈ 1 + μ · D_L ≈ 15 vs the predicted 1 + μ = 3.

**Status.**  No closed-form expression for the true P(L = l) at d > 0 is
derived in Stein & Werner (2025) or elsewhere in this project.  The
compound-Poisson prediction labeled "Markovian theory" in plots is therefore
valid at d = 0 and increasingly inaccurate as d → b.  Improving the prediction
would require deriving the true divisional distribution for the birth–death
process; using the empirically measured leaf depths in place of the Poisson
approximation in Eq. 24 is possible in principle but not a theoretical result.

### Mean divisional depth (Theorem 1 / Corollary, Stein & Werner 2025)

```
predict_averageGeneration(b, d, N):

  d = 0:  ⟨l⟩ = 2(γ + ln N − 1)
  d > 0:  ⟨l⟩ = 2b·t_N − 2,   t_N = (γ + ln(N(b−d)/b)) / (b−d)

  γ = 0.5772156649...  (Euler–Mascheroni constant)
```

---

## Deterministic case (fixed division time)

When every cell divides at a fixed time 1/b after birth (no stochasticity in
timing), the tree is a perfect binary tree.  All N leaves share the same
divisional depth

$$l^* = b \cdot t$$

and the divisional distribution is a **point mass**: E[D_{l*}] = N, E[D_l] = 0
for l ≠ l*.

Theorem 3 (Eq. 24) collapses to the single l = l* term:

$$\boxed{E[M_j] = N \cdot \text{Pois}(j;\; m \cdot l^*) = N \cdot \text{Pois}(j;\; m b t)}$$

**Key consequences:**

- **Mean mutations per cell:** E[M/N] = m · b · t
- **Distribution:** pure Poisson — no overdispersion from the tree structure
- **Index of dispersion:** Var[M] / E[M] = 1
- **Mutation-rate estimator** `var/mean − 1 = 0`, regardless of m — the estimator
  fails completely for deterministic trees because there is no divisional variance
  to exploit

The deterministic and Markovian cases bracket the biologically relevant regime:
real cells have partially regulated (non-Markovian) division times, placing the
true scMB overdispersion between these two extremes.

---

## Summary: index of dispersion by timing model

| Model | Divisional distribution | Index of dispersion (Var[M]/E[M]) |
|-------|------------------------|-----------------------------------|
| Deterministic (fixed 1/b) | δ(l, l*) | 1 |
| Non-Markovian, Gamma(k) | 0 < σ²_l < ⟨l⟩ | 1 < · < 1 + m |
| Markovian (exponential) | Pois(⟨l⟩) | 1 + m |

The Gamma(k) tree is more balanced than the exponential tree (σ²_l ≈ ⟨l⟩/k for
large k), so its index of dispersion ≈ 1 + m/k — smaller than the Markovian
value by a factor of k.
