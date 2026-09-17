# `2_subsampling/` — stage 1b: uniform `n`-cell subsamples of each full tree

For every raw shard under `data/raw/<scenario>/<model>/`, draws `n` of the tree's
`N` leaf cells uniformly without replacement and serializes the *induced* tree —
the sampled leaves plus every ancestor of a sampled leaf — under
`data/raw_subsampled/<scenario>/<model>/`. Because the observables measured on the
subsamples are the same ones the full-tree processing already computes,
`2_subsampling/` is the source stage for everything that reads a "sample-like"
copy-number or SFS view of the simulations.

The draw itself is `MutationLoadDynamics.sample_leaves` (v0.3.0+). This stage
supplies the sample-size table, the seed derivation, the atomic write, and the
per-scenario shard drivers (`X_helpers/subsampling.jl`).

## Invariants this stage depends on

**Prune, do not collapse.** Every division ancestral to a sampled cell is still
a node in the induced tree — no unary node is merged into its child. This is why
`mutations_per_cell` and `leaf_depths` on the subsampled tree return exactly the
same per-cell values as on the full tree; collapsing would turn `leaf_depths`
into a count of bifurcations that survived sampling, a property of the sample
rather than of the cell. The invariant lives upstream, in `sample_leaves`, and
is tested in `MutationLoadDynamics.jl`'s own `test/sampling.jl`.

**Co-indexing.** Index `i` of any array under `data/raw_subsampled/` and (later)
`data/processed_subsampled/` is the same simulation as index `i` under
`data/raw/` and `data/processed/` — both stages drop `nothing`-tree simulations
in the same source order. Full-vs-sample comparisons are therefore a plain array
zip with no join key. Do not reorder or filter either side.

**Sample sizes are a fixed table**, not computed. `sample_sizes(N_target)` in
`X_helpers/subsampling.jl`:

| `N_target` | sizes (largest-first) |
|---|---|
| 1 000 | `[100]` |
| 1 024 | `[102]` |
| 10 000 | `[1000, 100]` |
| 16 384 | `[1638, 164]` |

An unknown `N_target` errors. Adding a population size means adding a row — new
sizes never appear silently.

**Writes are atomic.** Both this stage and `3_processing/subsampled_trees/`
decide "already done" from `isfile`, so `serialize_atomic` writes to `<path>.tmp`
and renames. Without it, an interrupted run leaves a truncated-but-present file
that every subsequent run silently skips. Use it for any new output here.

**`SubsampleResult` has no draw-replicate field.** The design is one draw per
`(simulation, n)`. Supporting *k* independent draws per simulation would need
either a new struct field — breaking the 525 existing shards — or overloading
`n` in the filename. Decide that before adding it, not after.

## Scripts

| Script | Scenario | Source | Output |
|---|---|---|---|
| `subsample_neutral.jl` | `neutral` | `data/raw/neutral/` | `data/raw_subsampled/neutral/` |
| `subsample_sel1.jl` | `selection_1` | `data/raw/selection_1/` | `data/raw_subsampled/selection_1/` |
| `subsample_sel2.jl` | `selection_2` | `data/raw/selection_2/` | `data/raw_subsampled/selection_2/` |
| `verify_subsampling.jl` | — | — | unit suite (164 assertions) |
| `inventory_subsampled.jl` | — | — | filename audit |

Each `subsample_*.jl` mirrors the grid of the paired `1_simulation_runs/`
scripts exactly — same `N`, `d`, `s` loops, same interpolated filenames — and
delegates the actual work to `subsample_file(raw_path, out_dir,
sample_sizes(N))`.

Output filenames append `_n<n>`:

- `neutral_gamma_N10000_d0.5_k5.0_n100.jls` (from `_k5.0.jls`)
- `sel1_gamma_N10000_d0.5_k5.0_s0.3_n100.jls`
- `sel2_gamma_N10000_d0.5_k5.0_s0.15_M10.0_n100.jls`

## What has been run

- **525 shards, 4.2 GB on disk.**
- Neutral: 14 shards (6 gamma + 6 markov + 2 deterministic; one file per
  `N_target` at each entry of `sample_sizes(N)`).
- Selection 1: 280 raw shards × their sample sizes (gamma/markov `N ∈ {1000,
  10000}` gives one sample size each, deterministic `N ∈ {1024, 16384}` gives
  one each) = one subsampled shard per raw shard's `n`-list entry.
- Selection 2: 56 raw shards × ditto.

## Tests and audit

- `verify_subsampling.jl` — the unit suite. Run as
  `julia --project=. -t auto code/1_data_generation/2_subsampling/verify_subsampling.jl`.
  Its strongest test draws `n = N_full` and asserts the induced tree reproduces
  the arrays under `data/processed/` exactly.
- `inventory_subsampled.jl` — rebuilds every expected filename from the
  parameter grids and compares against disk, so a naming drift surfaces as a
  missing file rather than a silently-skipped shard. `--deep` additionally
  reads back every generated file (minutes, ~4.75 GB).

## Running

```bash
julia --project=. -t auto code/1_data_generation/2_subsampling/<script>.jl
```

Files sit at the top of `2_subsampling/`, one directory deeper than
`1_data_generation/` — so the `ROOT` header uses one fewer `dirname` than the
`1_simulation_runs/*/*.jl` scripts.
