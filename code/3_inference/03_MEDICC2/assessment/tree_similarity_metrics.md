# Tree similarity metrics — options for MEDICC2 vs. truth comparison

Notes drafted 2026-09-10 in response to "is there a systematic way to
investigate topology similarity between inferred and real trees?".

## Classical topology / branch-length metrics

- **Robinson–Foulds (RF)** — the workhorse. Counts bipartitions in one tree
  not present in the other. Normalized RF divides by max ($2n - 6$).
  Available as `ape::RF.dist()`. Fast, everyone knows it.
  **Weaknesses:** saturates quickly (one wrong split → jumps a lot), gives
  the same score for very different disagreements, ignores branch lengths.

- **Weighted RF** (Kuhner–Felsenstein) — RF with $\sum |\ell_i - \ell'_i|$
  over matched splits. `phangorn::wRF.dist()`.
  **This is exactly what a parsimony-saturation hypothesis asks for**: it
  separates "MEDICC2 got the topology roughly right but compressed the long
  branches" from "MEDICC2 got the topology wrong".

- **Quartet distance** — count 4-taxon subsets with disagreeing topology.
  Better distributed than RF, more informative near the max.
  `TreeDist::QuartetDivergence()`.

- **Cophenetic correlation** — Pearson (or Spearman) correlation of the two
  trees' pairwise leaf-to-leaf tree distances. Continuous, branch-length
  aware, intuitive. Compute with
  `cor(cophenetic.phylo(t1), cophenetic.phylo(t2))`.

- **Path difference distance** (Steel–Penny) — L2 norm of pairwise
  path-length differences.

- **Kendall–Colijn** — a vector representation of each tree from
  internal-node paths, then Euclidean distance.
  `TreeDist::KendallColijn()`. Popular in recent phylogenetics work
  because it interpolates smoothly between topology-only and length-only
  comparisons.

## Cell-relationship / "did we get who's related to whom right" metrics

These usually match your intuition better than RF for single-cell trees:

- **Ancestor–Descendant (AD) matrix agreement.** For every ordered pair
  $(i, j)$: is $i$ an ancestor of $j$ in the tree? Build the
  $n \times n$ boolean matrix from truth and MEDICC2; compare with
  precision / recall / F1 or Jaccard. This is *the* standard in the tumor
  single-cell phylogeny benchmarks (SCITE, SCARLET, SPhyR papers all
  report it).

- **Common Ancestor Set (CASet) distance** (DiNardo et al. 2019,
  *Bioinformatics*) — for each pair $(i, j)$, take the set of mutations on
  the common ancestor's edge; compare set overlap. Designed specifically
  for tumor phylogeny comparison.

- **Distinctly Inherited Set (DIS/DISC)** — companion metric to CASet.
  Focuses on mutations inherited by $i$ but not $j$.

- **MP3** (Ciccolella et al. 2021, *Bioinformatics*) — comparison via
  multi-labeled trees with mutation similarity weights. Handles the
  ancestral-CN case well.

- **Triplet distance** (rooted analogue of quartet) — for each 3-taxon
  subset, which of the 4 possible rooted topologies did each tree pick?
  `TreeDist::TripletDistance()`.

- **Sister precision/recall** — of the cherries (adjacent leaf pairs) in
  truth, what fraction survive in MEDICC2?

## What MEDICC2 benchmarks typically use

The MEDICC2 paper (Kaufmann et al. 2022, *Genome Biology*) evaluated
against simulated ground truth using **normalized RF** and by comparing
**pairwise MEDICC2 distances against true divisional distances**
(essentially the cophenetic correlation, in your language). They didn't
use CASet/MP3 because those need mutation-set semantics that CN events
don't have out of the box.

Comparable methods (SCARLET, HATCHet, ReMixT) mostly report
**AD accuracy + normalized RF + a branch-length or event-count
correlation**.

## Recommendation for this repo

For the saturation hypothesis, four quantities together cover topology and
length independently:

1. **Normalized RF** — one number, everyone recognizes it.
2. **Weighted RF** — how much of the tree "distance" is compression vs.
   real topology error.
3. **Cophenetic correlation** — Pearson $r$ between the two $\binom{n}{2}$
   leaf-pair distance vectors; a saturation-only signal shows up as $r$
   close to 1 in *rank* but a slope < 1 in the linear fit.
4. **Ancestor–Descendant F1** — the "did we get relationships right?"
   number, directly interpretable.

You already have both trees loaded, and `TreeDist` + `ape` cover all four
with a few lines. Natural fit: add these as a fourth panel to
`plot_saturation.R` (or a new `plot_topology_metrics.R` per the existing
"one script per figure" convention) so it lives next to the burden
diagnostics.

## Open question before implementation

If we go ahead, one design choice matters:

- **Per-simulation numbers only** (matches today's scope) — one small
  panel of four numbers per sim, next to the burden diagnostics.
- **Aggregate across sims** — a scatter of $(\text{E}[M], D,
  \text{RF}, r_\text{coph})$ per completed sim, the natural next step
  once more of the sweep completes.
