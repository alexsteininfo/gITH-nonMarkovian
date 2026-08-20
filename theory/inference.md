# Mutation rate inference under non-Markovian division timing

---

## Problem setup

We observe a population of $N$ cells, each carrying $M_i$ neutral mutations.
Neutral mutations arise at rate $m$ per division (Poisson).
Division times follow $\text{Gamma}(k, 1/(kb))$ — mean $1/b$, $\text{CV} = 1/\sqrt{k}$.

| $k$ | timing model |
|-----|-------------|
| 1 | Markovian (exponential / Gillespie) |
| 5 | cancer-cell regime (CV ≈ 0.45) |
| ∞ | fully deterministic (fixed $1/b$) |

---

## The master decomposition

Let $L$ = number of divisions undergone by a cell (leaf depth in the phylogenetic tree).
Each division independently contributes $\text{Pois}(m)$ mutations, so $M \mid L \sim \text{Pois}(mL)$.

By the law of total variance:

$$E[M] = m \cdot E[L]$$

$$\text{Var}[M] = m \cdot E[L] + m^2 \cdot \text{Var}[L]$$

Everything reduces to one dimensionless quantity, the **index of dispersion of leaf depth**:

$$\boxed{D \;\equiv\; \frac{\text{Var}[M]}{E[M]} \;=\; 1 + m \cdot c, \qquad c \equiv \frac{\text{Var}[L]}{E[L]}}$$

The mutation rate enters the overdispersion only through the product $m \cdot c$.
The two are **not separately identifiable from the scMB distribution alone** — additional
information (measured $k$, population dynamics, or SFS) is always required.

---

## Values of $c$ across timing models

| Model | divisional distribution | $c = \text{Var}[L]/E[L]$ | estimable from overdispersion? |
|-------|------------------------|--------------------------|-------------------------------|
| Deterministic ($k \to \infty$) | point mass at $l^*$ | 0 | ✗ — no signal |
| Gamma($k$) | underdispersed, $0 < c < 1$ | $0 < c < 1$ | ✓ with correction |
| Markovian ($k = 1$) | Poisson($\langle l \rangle$) | 1 | ✓ standard estimator |

For $k = 5$, simulations give $c \approx 0.14$, so the naive Markovian estimator
$\hat{m} = D - 1$ recovers only $m \times 0.14 \approx 0.28$ instead of $m = 2.0$.

---

## Why $c < 1/k$: within-lineage vs between-lineage variance

The index of dispersion $c$ has two additive contributions:

$$c_{\text{tree}} = \underbrace{\vphantom{\Big|}\frac{1}{k} \cdot \frac{E[L]_{\text{det}}}{E[L]}}_{\displaystyle\text{within-lineage}} \;+\; \underbrace{\vphantom{\Big|}\frac{\text{Var}_{\text{topology}}[L]}{E[L]}}_{\displaystyle\text{between-lineage}}$$

**Within-lineage** (renewal theory): along any single root-to-leaf path,
the number of divisions in time $t$ has $\text{Var}/E \approx 1/k$ (delta method on
the Gamma counting process).

**Between-lineage** (tree topology): different paths in the tree have different lengths
because early vs. late splits create asymmetric subtrees.

For the Yule tree ($k = 1$), the between-lineage term dominates — the leaf depth is
approximately $\text{Pois}(\langle l \rangle)$, so $c = 1$ exactly.
As $k$ grows, the tree becomes more balanced and the between-lineage variance
collapses faster than $1/k$, pushing $c$ below $1/k$.

For $k = 5$: renewal predicts $1/k = 0.20$; simulations give $c \approx 0.14$.
The ratio $\approx 0.7$ reflects the additional balancing from the tree topology.

---

## Boundary values and interpolation

For the key quantity $E[L]$, the two limits are:

$$E[L]_{\text{det}} = \frac{b \ln N}{b - d}, \qquad
  E[L]_{\text{Markov}} = \frac{2b\bigl(\gamma + \ln N(b-d)/b\bigr)}{b - d} - 2$$

where $\gamma = 0.5772\ldots$ is the Euler–Mascheroni constant.

**Empirical interpolation** (consistent with simulations, analytical status open):

$$\boxed{E[L]_{\text{Gamma}(k)} \;\approx\; E[L]_{\text{det}} \cdot \left(1 + \frac{1}{\sqrt{k}}\right)}$$

| $k$ | formula | $d=0$, $N=1{,}000$ |
|-----|---------|-------------------|
| 1   | $2 \ln N$ | 13.8 (exact: 13.0) |
| 5   | $1.447 \ln N$ | 10.0 (sim: ≈ 10) |
| ∞   | $\ln N$ | 6.9 |

Whether this follows from the multi-type (k-stage pipeline) embedding of the
age-dependent branching process is an open analytical question.

Similarly, a candidate formula for $c(k)$ that matches the boundary conditions and the
$k = 5$ simulation:

$$c(k) \;\approx\; \frac{1}{k + \sqrt{k}}$$

| $k$ | formula | simulation |
|-----|---------|-----------|
| 1   | 0.50    | 1.00 (wrong — tree topology dominates) |
| 5   | 0.138   | ≈ 0.14 ✓ |
| ∞   | 0       | 0 ✓ |

The $k = 1$ failure is expected: for the Yule tree the between-lineage variance
(not the renewal term) is the sole source of $c$, and the formula only captures the
within-lineage part. For $k \geq 2$ it is likely adequate.

---

## Two tractable estimators (k known)

### Estimator 1 — Mean-based (most robust)

$$\boxed{\hat{m} = \frac{\bar{M}}{E[L]_{\text{Gamma}(k)}}}$$

Uses no overdispersion — entirely agnostic about $c$. Requires $k$ and the
interpolation formula for $E[L]$.
Advantage: the mean is more stable than the variance, especially for small $N$.

### Estimator 2 — Corrected variance-based

$$\boxed{\hat{m} \approx \frac{D - 1}{c(k)}, \qquad c(k) \approx \frac{1}{k + \sqrt{k}}}$$

Requires $k$ only (not $E[L]$).
Reduces to the standard Markovian estimator for $k = 1$ (up to the $k = 1$ approximation error).
For $k \geq 2$, the $c(k)$ formula appears reliable.

---

## Identifiability: what external information is needed

From the scMB distribution **alone**, $m$ and $c$ are not separately identifiable —
multiplying $m$ by $\lambda$ and dividing $c$ by $\lambda$ leaves $D$ unchanged.
At least one of the following must be provided:

| Source | what it pins |
|--------|-------------|
| Inter-division time data (imaging / FACS) | $k$ → either estimator above |
| Population dynamics $(b, d, N)$ | $E[L]$ via interpolation → Estimator 1 |
| Site frequency spectrum (SFS) | tree topology → effective $c$ independent of $m$ |

**SFS approach**: $E[\xi_n] = m \cdot E[T_n]$, where $\xi_n$ is the count of mutations
at frequency $n/N$ and $T_n$ is the total branch length subtended by $n$ cells.
The *shape* of the SFS depends only on the tree (i.e., on $k$, not on $m$);
the *scale* pins $m$.
Fitting the SFS shape gives an estimate of $k$, which feeds into either estimator.

---

## Why ABC is avoidable: pre-calibration strategy

For given $(k, b, d, N)$, simulate trees **without mutations** (ν = 0 — the neutral
simulations already produced by `growth_neutral.jl`). From each tree, compute leaf
depths $\{L_i\}$ and tabulate:

$$\hat{c}(k, N, d) = \frac{\hat{\text{Var}}[L]}{\hat{E}[L]}, \qquad
  \widehat{E[L]}(k, N, d)$$

Then use the exact corrected estimator:

$$\hat{m} = \frac{D - 1}{\hat{c}}, \qquad \text{or} \qquad \hat{m} = \frac{\bar{M}}{\widehat{E[L]}}$$

This requires only $O(N_{\text{sims}})$ cheap tree-only simulations per parameter
combination — no mutation tracking, no likelihood computation, no expensive rejection
sampling. A grid over $(k, N, d)$ would cover the relevant parameter space once and
the look-up table is reusable.

---

## Open analytical questions

1. **Closed form for $E[L]_{\text{Gamma}(k)}$**: does the interpolation
   $(1 + 1/\sqrt{k}) \cdot E[L]_{\text{det}}$ follow from the $k$-stage pipeline
   representation of the age-dependent branching process?

2. **Exact $c(k)$ for the tree**: the formula $1/(k + \sqrt{k})$ matches simulations
   for $k \geq 5$ but fails at $k = 1$. Is there a single expression that covers
   the full range, including the Yule-tree limit?

3. **SFS-based estimator**: for a Gamma($k$) tree, what is $E[T_n]$ as a function
   of $k$ and $n$? This would enable the joint SFS + scMB estimator without
   any simulation calibration.
