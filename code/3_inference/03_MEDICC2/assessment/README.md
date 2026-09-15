# MEDICC2 assessment

Diagnostic figures comparing the ground-truth simulation output to the MEDICC2
reconstruction, per simulation. Implements the design at
`docs/superpowers/specs/2026-09-10-medicc2-assessment-design.md`.

## Scripts

Each script produces one PNG under
`figures/3_inference/03_MEDICC2/assessment/`. All take `--sim <relpath>`
(default: `neutral/gamma/neutral_gamma_N1000_d0.5_k5.0_n100/sim1`).

- `plot_cn_profiles.R` — allele-specific CN profiles for 6 cells (4 stratified
  by truth burden, 2 random) with truth overlaid on MEDICC2. **(task A)**
- `plot_phylogenies.R` — writes two PNGs per sim: `phylogenies_paired__*.png`
  (two ggtree panels with matched tip order) and `phylogenies_tanglegram__*.png`
  (ape::cophyloplot with connecting lines). **(task B)**
- `plot_saturation.R` — three-panel PNG: branch-length histograms (truth vs
  MEDICC2), per-cell burden histograms, moments table (`E[M]`, `Var[M]`,
  `D=Var/E`, `D-1`, with a MEDICC2/Truth ratio column). The `D-1` ratio is
  the headline diagnostic of parsimony saturation. **(task C)**
- `plot_topology_metrics.R` — tree-topology + branch-length agreement in one
  PNG: cophenetic pairwise-distance scatter (both trees in event units) with
  a fitted regression line, plus a summary table of ten metrics — normalised
  and raw Robinson-Foulds, weighted RF (raw + normalised), cophenetic
  Pearson / Spearman / slope, and sister-pair precision / recall / F1. Split
  from `plot_saturation.R` per the "one script per figure" convention; see
  `tree_similarity_metrics.md` for the rationale behind each metric.

## Helper

All disk I/O and tree walking lives in `code/X_helpers/medicc2_io.R`. Its
loaders take an absolute directory path (`truth_dir` = `data/CN_subsampled/...`,
`medicc2_dir` = `data/MEDICC2/treeinference/...`); callers do
their own path resolution.

Smoke test: `code/X_helpers/verify_medicc2_io.R` (run once after helper
edits).

## Running

```
micromamba run -n R Rscript code/3_inference/03_MEDICC2/assessment/plot_cn_profiles.R
micromamba run -n R Rscript code/3_inference/03_MEDICC2/assessment/plot_phylogenies.R
micromamba run -n R Rscript code/3_inference/03_MEDICC2/assessment/plot_saturation.R
micromamba run -n R Rscript code/3_inference/03_MEDICC2/assessment/plot_topology_metrics.R
```

To target a different sim:

```
... plot_cn_profiles.R --sim selection_1/gamma/sel1_gamma_N1000_d0.0_k5.0_s1.0_n100/sim1
```

## How this differs from the MEDICC2 paper's own validation

Reference: Kaufmann et al. 2022, *Genome Biol* 23:241,
[DOI: 10.1186/s13059-022-02794-9](https://doi.org/10.1186/s13059-022-02794-9)
(source: PubMed).

### What the MEDICC2 paper does

Verbatim from the paper's "Simulating genome evolution" methods section:

- Tree topology: "first a tree topology for a given number of leaves was
  created by **randomly joining sample labels** and rooting the tree at the
  diploid."
- Branch lengths: "The branch lengths and therefore the number of events per
  branch were determined using a **Poisson distribution** with λ = μ·l·r"
  (μ = 1, l = 440 segments, r ∈ {0.01, 0.025, 0.05}).
- Genome: 2 × 22 chromosomes × 10 uniform segments = 440 segments. Events
  include chromosomal / focal gains and losses, breakage-fusion-bridge, WGD,
  balanced and unbalanced translocations, inversions.
- Sample sizes: leaves ∈ [5, 10, 15, 20, 50, 100, 250, 500]; **leaves = the
  full population** (no subsampling).
- Metrics: **generalized Robinson-Foulds (GRF)**, regular Robinson-Foulds,
  Quartet distance.
- Result: "MEDICC2 outperforms other methods for all ranges of mutation rates
  and tree sizes, especially in the presence of WGDs".

### Head-to-head

| Aspect | MEDICC2 paper (Kaufmann 2022) | This repo |
|---|---|---|
| Tree topology source | Random uniform ("randomly joining sample labels") — a timeless PDA-like tree | Full lineage tree from a **birth-death simulation with actual division dynamics** (`MutationLoadDynamics.jl`) |
| Branch-length model | Poisson(μ·l·r) assigned directly to each branch | **Emergent**: # events per edge is a Poisson process along wall-clock time × edge duration |
| Division timing | Not modeled — trees are timeless with Poisson-length branches | **Deterministic / gamma(k=5) / exponential (Markov)** — three distinct waiting-time clocks |
| Death rate `d` | Not modeled | Explicit `d ∈ {0.0, 0.5, 0.9}` |
| Selection | Not modeled | Neutral + two selection scenarios (sel1: additive fixed effect; sel2: max-random) |
| Genome model | 440 segments; gains, losses, WGD, BFB, translocations, inversions | ~3000 segments (1 Mb bins); gains and losses only per the `sim1_truth_events.tsv` schema |
| Sample size | Leaves = full population, up to 500 leaves | **Uniform random subsample** of n=100 (or 102/164) from population N ∈ {1000, 1024, 10000, 16384} |
| Sampling depth | 100% | **10%, 1%, 0.6%** — much unobserved evolution collapsed onto sample-tree edges |
| Mutation rate range | μ ∈ [0.01, 0.025, 0.05] → mean events/branch ≈ 4.4 to 22 | Comparable scale by coincidence (target sim: max truth branch = 22 events, mean ≈ 5.8) |
| Evaluation metric | GRF + RF + Quartet on topology | Currently: branch-length + burden distributions + moments. Planned add-on: normalised RF, weighted RF, cophenetic-r, ancestor-descendant F1 (see `tree_similarity_metrics.md`) |
| Study of parsimony saturation | Not addressed as such | **The point of the work here** |
| Study of subsampling effect | Not addressed | Central design decision |
| Study of underlying division dynamics | Not addressed | Central design decision |

### What is (partly) duplicative

- Pure **topology-recovery accuracy** on Poisson-branch trees is what the
  MEDICC2 authors have already tested and where MEDICC2 wins.
- The general workflow "simulate → apply MEDICC2 → compare" is their setup.

### What is genuinely new

1. **Non-Markovian division timing.** MEDICC2's own simulator uses a
   memoryless Poisson process for branch lengths. The gamma(k=5) and
   deterministic clocks in this repo produce fundamentally different tree
   geometries (edge-length distributions, root-to-leaf variance, tree
   balance). MEDICC2 has never been tested on trees whose branch lengths
   are not drawn from an exponential/Poisson clock.
2. **Subsampling with hidden branches.** Sampling 100 leaves from
   1000–16384 collapses many hidden divisions onto each sample-tree edge.
   This is where parsimony saturation should be worst: a single "long"
   sample-tree edge may hide dozens of true CN events that MEDICC2
   compresses into a few parsimonious ones. MEDICC2 authors used the full
   population, so they never saw this regime.
3. **Death rate as an axis.** `d = 0.9` produces very unbalanced trees with
   very long unbroken branches — again the parsimony-worst-case. Nothing
   analogous in the paper.
4. **Burden-distribution moments, not just topology.** The paper reports
   RF-family accuracy on the topology only. This repo asks "did we get
   `E[M]` and `Var[M]` right?" — a prerequisite for the mutation-rate
   estimator in [`../../../../theory/inference.md`](../../../../theory/inference.md).
   Preliminary result on `neutral_gamma_N1000_d0.5_k5.0_n100/sim1`: MEDICC2
   preserves topology broadly but loses ~28% of total events
   (823 vs 1144) and ~46% of the per-cell burden variance
   (13 vs 24.1) — a specific, quantitative bias the paper's validation was
   structurally blind to.
5. **Selection scenarios.** Not covered in the paper.

### One-line pitch

MEDICC2's validation shows the tool recovers its own generative model; this
repo asks how well it holds up when the operating conditions match actual
tumor single-cell phylogenetics (birth-death dynamics, realistic clocks,
uniform-random subsampling of a small fraction of the population), with
downstream focus on burden-distribution moments rather than tree topology.
