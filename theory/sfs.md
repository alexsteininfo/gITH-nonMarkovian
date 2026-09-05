# Site Frequency Spectrum

## Definition

The site frequency spectrum (SFS) counts how many distinct mutations are carried by exactly $k$ cells in a final population of $N$ cells:

$$S_k = \#\{\text{mutations found in exactly } k \text{ cells}\}, \quad k = 1, \ldots, N.$$

Each mutation on a branch of the phylogenetic tree that subtends exactly $k$ leaves contributes 1 to $S_k$. With $\nu$ mutations accumulating per cell division (Poisson-distributed), the expected SFS equals the expected number of branches subtending $k$ leaves, scaled by $\nu$.

---

## Markovian birth–death model

**Reference:** Stein & Werner (2025), *Genetics* 230(4) iyaf101, Eq. 9 (originally Gunnarsson et al.).

For a population growing from 1 to $N$ cells under birth rate $b$ and death rate $d$, with $\nu$ neutral mutations per division, let $\omega = 2\nu$. The factor of 2 counts mutations on both internal branches (shared by $k \geq 2$ cells) and leaf branches (private, $k = 1$).

$$\mathbb{E}[S_k] = \omega N \sum_{j=0}^{\infty} \frac{(d/b)^j}{(k+j)(k+j+1)}$$

For $d = 0$ (pure birth / Yule process) the series collapses to:

$$\mathbb{E}[S_k] = \frac{\omega N}{k(k+1)}$$

**Implementation:** `predict_sfs_theory(k, b, d, N, mu)` in `analysis/helpers/theory.jl`.

### Total mutational burden

Summing over all frequencies:

$$\mathrm{TMB} = \sum_{k=1}^{N} \mathbb{E}[S_k] = \begin{cases} \omega(N-1) & d = 0 \\ \omega N \cdot \frac{-\ln(\sigma - \rho/N)}{\rho} & d > 0 \end{cases}$$

where $\rho = d/b$ and $\sigma = 1 - \rho$.

**Implementation:** `predict_tmb_theory(b, d, N, mu)` in `analysis/helpers/theory.jl`.

---

## Deterministic model (perfect binary tree)

When all cells divide synchronously (Dirac waiting times, division time $= 1/b$) and $N = 2^L$ is a power of two, the phylogeny is a **perfect binary tree** of depth $L = \log_2 N$.

At depth $d$ from the root there are $2^d$ branches, each subtending exactly $N / 2^d$ leaves. With $\nu$ mutations per division:

$$\mathbb{E}[S_k] = \frac{N \nu}{k}, \quad k \in \{1,\, 2,\, 4,\, \ldots,\, N/2\}$$

$$\mathbb{E}[S_k] = 0 \quad \text{otherwise.}$$

The SFS is nonzero **only at power-of-2 frequencies**. This is the direct fingerprint of synchronous cell division: all branches at a given depth subtend the same number of leaves.

**Implementation:** `predict_sfs_deterministic(k, N, nu)` in `analysis/helpers/theory.jl`.

### Total mutational burden

$$\mathrm{TMB} = \sum_{j=0}^{L-1} \frac{N\nu}{2^j} = N\nu \cdot \frac{2(1 - N^{-1})}{1} = 2(N-1)\nu$$

This equals the Markovian TMB at $d = 0$ ($\omega(N-1) = 2\nu(N-1)$), confirming that total mutations are conserved — only the distribution across frequencies differs.

### Comparison with the Markovian prediction

The Markovian formula (Eq. 9) distributes the budget $2(N-1)\nu$ across all $N-1$ frequency bins with $1/[k(k+1)]$ weighting, giving:

$$\mathbb{E}^{\mathrm{Markov}}[S_k] = \frac{\omega N}{k(k+1)} \approx \frac{\omega N}{k^2} \quad (k \ll N)$$

The deterministic model concentrates the same budget into only $L = \log_2 N$ bins. At a power-of-2 frequency $k$, the ratio is:

$$\frac{\mathbb{E}^{\mathrm{det}}[S_k]}{\mathbb{E}^{\mathrm{Markov}}[S_k]} \approx \frac{N\nu/k}{\omega N / k^2} = \frac{k}{2} \quad (k \ll N)$$

For small $k$ (relevant part of the spectrum), the ratio grows linearly with $k$ — i.e., the deterministic SFS is not simply shifted up by a constant factor. At $k=1$ the Markovian formula gives $\omega N / 2$ and the deterministic gives $N\nu = \omega N / 2$, so they **agree exactly at $k = 1$**. The discrepancy grows for larger $k$: the Markovian formula decays as $1/k^2$ while the deterministic decays as $1/k$, so the deterministic SFS is increasingly higher at intermediate frequencies.

---

## Subsampling: the hypergeometric projection

Every formula above describes the SFS of a **complete** population of $N$ cells.
Real single-cell data never is: one sequences $n \ll N$ cells. This section gives
the exact map from the full-population SFS to the SFS of a uniform sample, and two
consequences that matter for interpreting sampled data.

### The projection

Sample $n$ of the $N$ cells uniformly without replacement. A mutation carried by
exactly $j$ of the $N$ cells is carried by exactly $k$ of the $n$ sampled cells with
the hypergeometric probability

$$P(k \mid j) = \frac{\binom{j}{k}\binom{N-j}{n-k}}{\binom{N}{n}},$$

since sampling $n$ cells from a population in which $j$ are carriers is exactly a
hypergeometric draw. Summing over the mutations at each full-population frequency
gives the expected sample SFS:

$$\boxed{\;\mathbb{E}\!\left[S_k^{(n)}\right] \;=\; \sum_{j=k}^{N-n+k} \mathbb{E}\!\left[S_j^{(N)}\right] \frac{\binom{j}{k}\binom{N-j}{n-k}}{\binom{N}{n}}, \qquad k = 1, \ldots, n.\;}$$

Three properties are worth stating explicitly.

**It is linear and model-independent.** The projection is a fixed matrix
$H_{kj} = P(k \mid j)$ depending only on $(N, n)$ — never on the division-time
distribution. So the *same* operator applies to the Markovian, gamma and
deterministic spectra, and any difference between the sampled spectra is inherited
entirely from the difference between the full spectra. This is what makes the
projection usable as a correctness oracle: it predicts the sampled SFS from the full
SFS with no free parameters and no reference to the timing model.

**Mutations are lost, and the loss is frequency-dependent.** A mutation lands in the
unobserved bin $k = 0$ with probability $\binom{N-j}{n}\big/\binom{N}{n}$. For a
singleton ($j = 1$) this is $1 - n/N$: only the fraction $n/N$ of private mutations
survive sampling at all. High-frequency mutations are almost never lost. The total
number of segregating sites therefore falls, and it falls hardest at the
low-frequency end where most of the sites are.

**Per-cell burden is conserved exactly.** Because the induced tree retains every
ancestor of every sampled cell (see `analysis/helpers/subsampling.jl`), a sampled
cell carries exactly the mutations it carried in the full tree. Equivalently, using
$\mathbb{E}[k \mid j] = nj/N$,

$$\sum_{k=1}^{n} k \,\mathbb{E}\!\left[S_k^{(n)}\right] = \frac{n}{N} \sum_{j=1}^{N} j \,\mathbb{E}\!\left[S_j^{(N)}\right],$$

so the mean burden per cell, $\frac{1}{n}\sum_k k\,S_k^{(n)}$, is unchanged. The
single-cell mutational burden distribution is thus *invariant* under sampling — only
the SFS is distorted. That asymmetry is why `mut_per_cell` and `leaf_depths` serve
as the correctness check on the subsampling code while `sfs` is the observable of
interest.

### Consequence 1: the sample SFS is **not** Eq. 9 with $N \to n$

The tempting shortcut — fit the closed-form prediction with the population size
replaced by the sample size — is wrong, and not by a little. Projecting the
Markovian $d = 0$ spectrum $\mathbb{E}[S_j] = \omega N / [j(j+1)]$ and comparing
against $\omega n / [k(k+1)]$:

| $k$ | ratio at $n/N = 0.1$ | ratio at $n/N = 0.01$ |
|---|---|---|
| 1 | 3.47 | 7.39 |
| 2 | 2.18 | 2.81 |
| 3 | 1.74 | 1.96 |
| 5 | 1.41 | 1.49 |
| 10 | 1.19 | 1.22 |

The true sample SFS is **elevated** at low frequencies relative to the naive
substitution, because the full spectrum's many high-$j$ mutations get only partially
captured and pile up in the low-$k$ bins. The ratio depends on the sampling
*fraction* $n/N$ rather than on $N$ and $n$ separately (the $1/[j(j+1)]$ spectrum is
scale-free), and it decays toward 1 as $k$ grows.

**Practical consequence for inference:** fitting Eq. 9 to subsampled data with $N$
set to the sample size will bias the inferred parameters, most severely through the
singleton bin. Either fit the projected prediction, or restrict the fit to
intermediate $k$ where the ratio is near 1. See `theory/inference.md`.

### Consequence 2: sampling erases the deterministic fingerprint at low frequency

The deterministic model's signature is that the SFS is nonzero **only** at
power-of-2 frequencies. Sampling smears each spike into a hypergeometric bump: the
spike at $j$ lands at mean position

$$\bar{k} = \frac{nj}{N}, \qquad \mathrm{sd} = \sqrt{n\,\frac{j}{N}\left(1 - \frac{j}{N}\right)\frac{N-n}{N-1}} \;\approx\; \sqrt{\bar{k}} \quad (j \ll N,\ n \ll N).$$

Adjacent spikes at $j$ and $2j$ project to $\bar{k}$ and $2\bar{k}$, a separation of
$\bar{k}$, while each has width $\approx\sqrt{\bar{k}}$. The separation-to-width
ratio is therefore $\sqrt{\bar{k}}$, and two neighbouring spikes remain
distinguishable only when

$$\bar{k} = \frac{nj}{N} \gtrsim 4, \qquad \text{i.e.} \qquad j \gtrsim \frac{4N}{n}.$$

Since the spikes sit at $j = 1, 2, 4, \ldots, N/2$, the number that survive is

$$\#\{\text{resolvable spikes}\} \;\approx\; \log_2\!\left(\frac{n}{4}\right),$$

which depends on the **sample size alone** — the population size $N$ drops out
entirely. Growing the tumour does not help; only sequencing more cells does.

Numerically (spike resolved when mean/sd $> 2$):

| $N$ | $n$ | spikes present | smeared | resolved | $\lfloor\log_2(n/4)\rfloor$ |
|---|---|---|---|---|---|
| 1024 | 102 | $j = 1 \ldots 512$ (10) | $j \le 32$ | $j \ge 64$ (4) | 4 |
| 16384 | 1638 | $j = 1 \ldots 8192$ (14) | $j \le 32$ | $j \ge 64$ (8) | 8 |
| 16384 | 164 | $j = 1 \ldots 8192$ (14) | $j \le 256$ | $j \ge 512$ (5) | 5 |

So at 10% sampling of a 1024-cell tree, six of the ten spikes are gone and the
low-frequency end looks like a smooth continuum — precisely the region where the
Markovian and deterministic predictions are most often compared. The
synchronous-division fingerprint is not destroyed by sampling, but it retreats to
frequencies $k \gtrsim 4$, and a plot restricted to the first few bins will not show
it.

### Numerical verification

The projection was checked against the simulated data by projecting the *measured*
mean full-tree SFS from `data/processed/` and comparing with the *measured* mean
subsampled SFS from `data/processed_subsampled/` — no fitting, no free parameters:

| shard | $n$ | predicted total | observed total | rel. error |
|---|---|---|---|---|
| `neutral_gamma_N1000_d0.5_k5.0` | 100 | 1090.75 | 1089.71 | 0.096% |
| `neutral_markov_N10000_d0.9` | 100 | 10222.66 | 10210.23 | 0.122% |
| `neutral_deterministic_N1024` | 102 | 943.75 | 943.98 | 0.025% |
| `neutral_deterministic_N16384` | 1638 | 15210.43 | 15210.36 | 0.000% |

Per-bin agreement at $k = 1, 2, 3, 10$ is within Monte-Carlo noise (0.01–3%, largest
where the bin count is smallest). Because the projection assumes only a uniform
draw and correct frequency bookkeeping, this simultaneously confirms that the
subsampling draw is uniform without replacement, that $S_k^{(n)}$ genuinely counts
mutations carried by exactly $k$ of the $n$ sampled cells, and that the induced tree
retains exactly the right mutations.

**Implementation:** the subsampling pipeline is `analysis/subsampling/` (stage 1b) and
`analysis/processing_subsampled/` (stage 2b); sample sizes are the fixed table in
`sample_sizes` in `analysis/helpers/subsampling.jl`. The projection itself is not yet
a script in the repo — see `TODO.md`.
