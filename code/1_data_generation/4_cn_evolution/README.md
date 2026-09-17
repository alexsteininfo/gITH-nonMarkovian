# `4_cn_evolution/` — stage 4: allele-specific copy-number evolution on the subsampled trees

Simulates allele-specific copy-number alterations (CNAs) on the stage-1b
subsampled trees with **CopyNumberEvolution.jl** (a separate, unregistered,
actively-developed package by the same author, `Pkg.develop`ed from
`../CopyNumberEvolution.jl`), and writes **MEDICC2-ready** input plus ground
truth to `data/CN_subsampled/`. Design rationale and every scoping decision:
`docs/superpowers/specs/2026-09-07-cn-evolution-design.md`.

Unlike the earlier stages, this one does not read `data/processed/` or
`data/processed_subsampled/` at all — the only input is the `tree_root` inside a
`SubsampleResult`, which `data/raw_subsampled/` already carries.

## The copy-number evolution model

```julia
CNAModel(rate = FromEdgeMutations())
```

is passed to `simulate_cnas(tree, hg38(:female), CN_MODEL; seed = seed)`. All
other `CNAModel` components are left at the **package default**, and the
defaults already implement every constraint set for this first pass (focal only,
no WGD, clean diploid root). The table below lists every component and the exact
value used.

| Component | Setting | What it means |
|---|---|---|
| **Assembly** | `hg38(:female)` | GRCh38 human reference. Chromosome lengths from the UCSC `hg38.chrom.sizes` table (retrieved 2026-09-04). **`:female`** → **22 autosomes each at ploidy 2 + X at ploidy 2 + Y at ploidy 0** (total 46 haplotype slots). Switching to `hg38(:male)` would give 22 autosomes × 2 + X × 1 + Y × 1 (also 46 slots) and change which chromosomes are eligible for the target draw. |
| **Bin grid** | `BinGrid(hg38(:female), 1_000_000)` | 1 Mb fixed-width bins over the whole assembly (all autosomes, X, Y). Coarsened from the package's 500 kb DLP+ default to halve `cells.tsv` size at the same cell count for this first pass; ≈ 2 875 autosomal bins. 500 kb is the realistic target once the first pass has been inspected. |
| **`rate`** | `FromEdgeMutations()` | **Explicitly requested.** `p = 1.0` (its default): one recorded mutation becomes exactly one CNA, no extra draw. Requires `edge_mutations`, which the MutationLoadDynamics.jl bridge (`CopyNumberEvolutionMutationLoadDynamicsExt`) populates from `cell.mutations` on every edge. The root's `edge_mutations` is left `nothing` — the founder's own mutations are not translated. |
| **`target`** | `UniformChromosome()` | Package default. Uniform over eligible chromosomes, then uniform over that chromosome's haplotypes. **Uniform over chromosomes**, not over base pairs — chromosome 1 and chromosome 22 are equally likely, so short chromosomes see more events per base than long ones. |
| **`extent`** | `ExtentMixture()` — `p_chromosome = 0`, `p_arm = 0`, `lengthdist = LogUniform(1e5, 1e8)` | Package default, and **is** the "focal-only" requirement (no arm- or chromosome-level events). Focal event length is drawn log-uniformly from **100 kb – 100 Mb**, then a uniform start on the target chromosome, then truncated at the chromosome end (matching MEDICC2). |
| **`kind`** | `GainLoss(0.5)` | Package default. `p_gain = 0.5`, `delta = 1` — half gains, half losses, each of magnitude 1 copy. A loss to CN 0 *is* loss of heterozygosity; no separate event type is needed. |
| **`wgd`** | `NoWGD()` | Package default, and **is** the explicit "no WGD" requirement. |
| **`viability`** | `RejectAndRedraw(min_total_cn = 1, max_attempts = 100)` | Package default. Rejects any event that would drive the total copy number (summed across haplotypes) below 1 at any touched position — i.e. forbids homozygous deletion of any size. A conservative first-pass choice for real-data-shaped output; set `min_total_cn = 0` to allow biallelic loss. |
| **`initial`** | `Diploid()` | Package default. Root profile is diploid with no truncal alterations. Chosen over `TruncalCNAs(founder_mutations(root))` because the founder's own mutation count is not attributable to any edge and would need its own CNA-generation choice. |

Tree conversion goes through the package extension:
`CopyNumberEvolution.PhyloTree(sim.tree_root)`, active once both
`CopyNumberEvolution` and `MutationLoadDynamics` are loaded. This assigns
`edge_divisions = 1` and `edge_mutations = cell.mutations` per edge — exactly
what `FromEdgeMutations` needs.

Seeding is `hash((stem, sim_index, :cn_evolution))` masked to `Int64` (see
`X_helpers/cn_evolution.jl`), a pure function of values already on disk, so any
run reproduces exactly.

## Scripts and scope

| Script | Scenario | Timing models in scope | Shards |
|---|---|---|---|
| `cn_neutral.jl` | `neutral` | deterministic, gamma, markov | 14 |
| `cn_sel1.jl` | `selection_1` | **gamma only** | 120 |
| `cn_sel2.jl` | `selection_2` | **gamma only** | 24 |

**158 shards total.** Two scoping decisions bring the first pass down from a
potentially several-hundred-GB sweep:

1. **Only the smallest sample size of each shard family**, i.e.
   `last(sample_sizes(N))` — the existing `sample_sizes` table returns each
   `N_target`'s sizes largest-first, so this is `n = 100` (from `N = 1000` or
   `N = 10000`), `n = 102` (from `N = 1024`), or `n = 164` (from `N = 16384`).
   `n = 1000` and `n = 1638` shards are out of scope for this pass — dropping
   them saves ~10× disk, and 100–200 cells is the scale of a real single-cell
   WGS copy-number study.
2. **Only the `gamma` timing model for `selection_1` and `selection_2`.**
   `neutral` keeps all three models — it is by far the smallest scenario by
   shard count. Adding `markov`/`deterministic` to the two selection scripts is
   a copy of the existing gamma loop; deferred because tripling that cost
   before the first pass has been inspected wasn't worth it.

`selection_1`/`selection_2`'s `gamma` grid only uses `N ∈ {1000, 10000}`, so
all 144 selection shards are at `n = 100`; the `n = 102`/`164` shards occur only
in `neutral`'s deterministic model, one shard each.

Per-shard grids (as coded in each `cn_*.jl`):

- `cn_neutral.jl` — gamma/markov: `N ∈ {1000, 10000} × d ∈ {0, 0.5, 0.9}`;
  deterministic: `N ∈ {1024, 16384}`.
- `cn_sel1.jl` — gamma only: `N ∈ {1000, 10000} × d ∈ {0, 0.5, 0.9} × s ∈
  0.1:0.1:2.0` (20 driver strengths).
- `cn_sel2.jl` — gamma only: `N ∈ {1000, 10000} × d ∈ {0, 0.5, 0.9} × s ∈
  {0.05, 0.10, 0.15, 0.20}`, `M = 10.0`.

## Per-simulation scope and resumability

**`SIMS_PER_SHARD = 1`** (`X_helpers/cn_evolution.jl`) — one simulation per
shard, in stored order (the same order `raw_subsampled/` already filters and
preserves from the raw simulations, not a new random pick). Resumable per
`sim_index`: raising `SIMS_PER_SHARD` later only adds the newly-included indices
without re-touching what is already on disk. A simulation is treated as done
when all four of its output files exist.

No `serialize_atomic` — the four `write_*` calls are ordinary buffered writes
that either complete or throw. A killed run would leave a `cells.tsv` with the
wrong row count rather than something a plain `isfile` check would treat as
done.

## Output layout

```
data/CN_subsampled/<scenario>/<model>/<stem>/sim<sim_index>_cells.tsv
data/CN_subsampled/<scenario>/<model>/<stem>/sim<sim_index>_truth_profiles.tsv
data/CN_subsampled/<scenario>/<model>/<stem>/sim<sim_index>_truth_events.tsv
data/CN_subsampled/<scenario>/<model>/<stem>/sim<sim_index>_tree.nwk
```

Example:
`data/CN_subsampled/selection_1/gamma/sel1_gamma_N10000_d0.5_k5.0_s0.3_n100/sim7_cells.tsv`

| File | Content | Writer |
|---|---|---|
| `cells.tsv` | Leaves only, projected onto the 1 Mb bin grid — **the actual MEDICC2 input** | `write_medicc2` |
| `truth_profiles.tsv` | Every internal + leaf node's profile — ground truth (segment-based, not binned) | `write_profiles` |
| `truth_events.tsv` | Every CNA event with order and edge | `write_events` |
| `tree.nwk` | Topology, `branchlength = :divisions` | `write_newick` |

The truth files are cheap (segment-based, not projected onto the bin grid) and
are the documented point of later comparison against MEDICC2's own inferred
tree and branch lengths (see `CopyNumberEvolution.jl`'s `docs/src/interop.md`).

## Disk footprint

**2.1 GB measured** — 1.6 GB of that is `cells.tsv`; the truth files and
newick trees add negligibly. At `SIMS_PER_SHARD = 1` and 1 Mb bins,
`cells.tsv` rows ≈ `(n_cells + 1) × n_bins` with `n_bins ≈ 2 875`
(hg38 autosomes only, half the 500 kb-bin count).

## Verification

No `verify_*.jl` / `inventory_*.jl` yet — premature before `SIMS_PER_SHARD` and
the 1 Mb resolution are validated against real usage. `medicc2` itself was not
run against the output on the machine this stage was built on (not installed);
the format was checked structurally instead (header line, autosome-only chrom
column, copy numbers within MEDICC2's 0–8 range).

## Extending the first pass

Every scoping decision is a one-spot edit:

- `SIMS_PER_SHARD` — the constant in `X_helpers/cn_evolution.jl`.
- Larger sample sizes — drop `last(...)` in each `cn_*.jl` and iterate the
  whole `sample_sizes(N)` vector, the shape `subsample_*.jl` already shows.
- `markov` / `deterministic` for the selection scenarios — add the missing
  model loops to `cn_sel1.jl` / `cn_sel2.jl`, mirroring `cn_neutral.jl`.

Resumability means any of these extensions only adds the newly-included work
without re-touching what is on disk.

## Running

```bash
julia --project=. code/1_data_generation/4_cn_evolution/cn_neutral.jl
julia --project=. code/1_data_generation/4_cn_evolution/cn_sel1.jl
julia --project=. code/1_data_generation/4_cn_evolution/cn_sel2.jl
```

Single-threaded — `-t auto` is harmless but unnecessary. Files sit at the top of
`4_cn_evolution/`, one directory deeper than `1_data_generation/`, so the `ROOT`
header uses one fewer `dirname` than the `1_simulation_runs/*/*.jl` scripts.
