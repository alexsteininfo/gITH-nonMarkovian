# `3_processing/subsampled_trees/` — stage 2b: observables from the subsamples

The same quantities `full_trees/` extracts, recomputed on the subsampled trees
under `data/raw_subsampled/`, into `data/processed_subsampled/<scenario>/<quantity>/`.
`sfs` is the quantity sampling genuinely distorts (see `theory/sfs.md`'s
hypergeometric projection); `mut_per_cell` and `leaf_depths` are unchanged from
the full-tree values for the same cells, by construction of stage 1b (every
ancestor of a sampled cell is retained, so root-to-leaf paths are unchanged).

## Scripts and grids

| Script | Source | Extra quantities | Output |
|---|---|---|---|
| `process_neutral_subsampled.jl` | `data/raw_subsampled/neutral/{det,gamma,markov}/` | — | `data/processed_subsampled/neutral/…` |
| `process_sel1_subsampled.jl` | `data/raw_subsampled/selection_1/…` | `injection` (**copied from `data/processed/selection_1/injection/`**) | `data/processed_subsampled/selection_1/…` |
| `process_sel2_subsampled.jl` | `data/raw_subsampled/selection_2/…` | `leaf_fitness`, `n_restarts` (**latter copied from `data/processed/selection_2/n_restarts/`**) | `data/processed_subsampled/selection_2/…` |
| `check_neutral_subsampled.jl` | — | — | Cross-checks subsampled vs. full-tree arrays. No sel1/sel2 equivalent yet. |

Filenames mirror the subsampled shards (stem includes `_n<n>`), e.g.
`neutral_gamma_N10000_d0.5_k5.0_n100.jls`.

## What is in each quantity file

For a subsampled shard of size `n` with `n_sims` entries:

| Quantity | Type | Length | Notes |
|---|---|---|---|
| `params` | `Vector{<Scenario>Params}` | `n_sims` | inherited from the source simulation |
| `mut_per_cell` | `Vector{Vector{Int}}` | `n_sims` (each = `n`) | identical to the corresponding cells' full-tree values |
| `sfs` | `Vector{Vector{Int}}` | `n_sims` (each = `n`) | `sfs[k]` = mutations shared by exactly `k` of the `n` sampled cells |
| `leaf_depths` | `Vector{Vector{Int}}` | `n_sims` (each = `n`) | identical to the corresponding cells' full-tree values |
| `leaf_fitness` (sel2) | `Vector{Vector{Float64}}` | `n_sims` | co-indexed with `mut_per_cell` |
| `injection` (sel1) | `Vector{Sel1Injection}` | `n_sims` | **copied** from `data/processed/selection_1/injection/` — a hard error if that stage hasn't run |
| `n_restarts` (sel2) | `Vector{Int}` | `n_sims` | **copied** from `data/processed/selection_2/n_restarts/` — a hard error if that stage hasn't run |

`sfs` has length `n`, not `N` — the deepest bin `sfs[n]` holds the mutations
carried by every sampled cell, including the founder's own.

Each shard is atomically written with `serialize_atomic` (see
`X_helpers/subsampling.jl`), so the `isfile`-based resumability check cannot see
a truncated file left by a killed run.

## Dependency chain

**`process_sel{1,2}_subsampled.jl` require the corresponding
`full_trees/process_sel{1,2}.jl` to have run first.** They copy the per-simulation
`injection` (sel1) or `n_restarts` (sel2) scalars across from `data/processed/`
rather than re-deriving them from the raw trees — re-deriving would mean a second
full pass over several GB of trees for a handful of numbers per simulation. Both
scripts hard-error, naming the required script, if those inputs are missing or
misaligned.

## Co-indexing invariant

Index `i` of any array here is the same simulation as index `i` of the matching
array under `data/processed/`. Both stages drop `nothing`-tree simulations in the
same source order. Paired full-vs-sample comparison is a plain array zip — no
join key. Do not reorder or filter either side.

## Cross-check

`check_neutral_subsampled.jl` asserts the subsampled arrays against the
full-tree arrays they are co-indexed with. Values that must match cell for cell
(`mut_per_cell`, `leaf_depths` at the sampled indices) are asserted equal;
values that are sampled (`sfs`) are checked for length and legality. There is no
`sel1`/`sel2` equivalent yet — see `TODO.md`.

## Running

```bash
julia --project=. -t auto code/1_data_generation/3_processing/subsampled_trees/process_neutral_subsampled.jl
julia --project=. -t auto code/1_data_generation/3_processing/subsampled_trees/process_sel1_subsampled.jl
julia --project=. -t auto code/1_data_generation/3_processing/subsampled_trees/process_sel2_subsampled.jl
julia --project=. -t auto code/1_data_generation/3_processing/subsampled_trees/check_neutral_subsampled.jl
```
