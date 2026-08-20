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
