# `3_processing/full_trees/` — stage 2: observables from the full lineage trees

One script per scenario. Each deserializes the stage-1 raw trees for that
scenario and walks them (via **MutationLoadDynamics.jl**'s `mutations_per_cell`,
`sitefrequencyspectrum`, `leaf_depths`, and — for `selection_2` —
`leaf_fitness`) to extract the arrays every downstream stage consumes.

Output layout: **one quantity per subdirectory**, one `.jls` per source shard,
same stem as the raw shard:

```
data/processed/<scenario>/<quantity>/<same-stem>.jls
```

## Scripts, quantities, and grids

| Script | Source | Quantities | Shards processed |
|---|---|---|---|
| `process_neutral.jl` | `data/raw/neutral/{det,gamma,markov}/` | `params`, `mut_per_cell`, `sfs`, `leaf_depths` | 14 (6 gamma + 6 markov + 2 det) |
| `process_sel1.jl` | `data/raw/selection_1/{det,gamma,markov}/` | + `injection` | 280 (120 gamma + 120 markov + 40 det) |
| `process_sel2.jl` | `data/raw/selection_2/{det,gamma,markov}/` | + `leaf_fitness`, `n_restarts` | 56 (24 gamma + 24 markov + 8 det) |

Grids mirror `1_simulation_runs/*/*.jl` exactly (same interpolated filenames).
`process_neutral.jl` reads only `N1024` and `N16384` from
`data/raw/neutral/deterministic/`; the leftover `N1000`/`N10000` files there are
orphans from an earlier parameterization.

## What is in each quantity file

For a shard with `n_sims` simulations that all have a non-`nothing` `tree_root`:

| Quantity | Type | Length | Notes |
|---|---|---|---|
| `params` | `Vector{<Scenario>Params}` | `n_sims` | full parameter record per sim |
| `mut_per_cell` | `Vector{Vector{Int}}` | `n_sims` (each `= actual_N`) | one count per alive leaf |
| `sfs` | `Vector{Vector{Int}}` | `n_sims` (each `= actual_N`) | `sfs[k]` = mutations shared by exactly `k` cells |
| `leaf_depths` | `Vector{Vector{Int}}` | `n_sims` (each `= actual_N`) | divisions from founder to each leaf |
| `injection` (sel1) | `Vector{Sel1Injection}` | `n_sims` | `t_inject`, `N_at_inject`, `driver_cell_id`, `driver_clone_size`, `n_attempts`, `seed` |
| `leaf_fitness` (sel2) | `Vector{Vector{Float64}}` | `n_sims` | co-indexed with `mut_per_cell` (both walk `Leaves` in the same order) |
| `n_restarts` (sel2) | `Vector{Int}` | `n_sims` | extinctions before the accepted attempt |

`actual_N ≥ N_target` because every model stops at the first `popsize ≥
N_target`, and deterministic timing overshoots to the next power of two.

`sfs` is length `actual_N`, so `sfs[actual_N]` is the count of mutations
present in every alive cell — including any founder mutation.

Simulations with `tree_root === nothing` are dropped in source order and a warning
is printed. `subsample_*.jl` and the subsampled processors drop in the same order,
so index `i` here co-indexes with index `i` there.

## Resumability

- **`process_sel1.jl` and `process_sel2.jl` are resumable.** A shard whose every
  quantity file (`params`, `mut_per_cell`, `sfs`, `leaf_depths`, plus scenario
  extras) already exists is skipped. Delete a processed file to force a redo.
- **`process_neutral.jl` is not resumable.** `process_file` serializes
  unconditionally, so re-running rewrites all 433 MB of
  `data/processed/neutral/`. Never run it to "check" something.

## Running

```bash
julia --project=. -t auto code/1_data_generation/3_processing/full_trees/process_neutral.jl
julia --project=. -t auto code/1_data_generation/3_processing/full_trees/process_sel1.jl
julia --project=. -t auto code/1_data_generation/3_processing/full_trees/process_sel2.jl
```
