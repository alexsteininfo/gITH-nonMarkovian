# `3_processing/` — stage 2 / 2b: extract observable arrays from the trees

Walks the stage-1 raw trees (`data/raw/`) and the stage-1b subsampled trees
(`data/raw_subsampled/`), computes the observables everything downstream reads,
and serializes one array per quantity per shard. This is where SFS, single-cell
mutational burden, leaf depths, and (for the selection scenarios) fitness /
injection / restart metadata become plain `Vector`s that plotting and inference
code can `deserialize` in milliseconds without ever touching the trees.

The tree statistics (`mutations_per_cell`, `sitefrequencyspectrum`,
`leaf_depths`, `leaf_fitness`) all live in **MutationLoadDynamics.jl** (v0.3.0+);
this stage is thin glue that opens shards, calls the package, and lays the
results out on disk.

## Subfolders

| Subfolder | What it does |
|---|---|
| `full_trees/` | Extracts arrays from the full lineage trees under `data/raw/` into `data/processed/<scenario>/<quantity>/`. Three scripts, one per scenario. |
| `subsampled_trees/` | Extracts the same arrays (with `sfs` sampled) from the induced trees under `data/raw_subsampled/` into `data/processed_subsampled/<scenario>/<quantity>/`. Three scripts, plus a neutral cross-check. |
| `checks/` | Gate that verifies this repo's stored arrays and subsample draws still agree with what `MutationLoadDynamics.jl` computes today. |

## What is written where

For each shard, one `.jls` per quantity, stem shared with the source shard:

```
data/processed/<scenario>/<quantity>/<same-stem>.jls
data/processed_subsampled/<scenario>/<quantity>/<same-stem>_n<n>.jls
```

Quantities are scenario-specific:

| Scenario | Quantities |
|---|---|
| `neutral` | `params`, `mut_per_cell`, `sfs`, `leaf_depths` |
| `selection_1` | `params`, `mut_per_cell`, `sfs`, `leaf_depths`, `injection` |
| `selection_2` | `params`, `mut_per_cell`, `sfs`, `leaf_depths`, `leaf_fitness`, `n_restarts` |

`injection` (a `Sel1Injection`) and `n_restarts` (an `Int`) are per-simulation
scalars that would otherwise force a second full pass over multi-GB tree files
to recover a handful of numbers per shard. They are stored here so the plotting
stage never re-opens the raw trees.

## Data volumes

| Directory | Size | Content |
|---|---|---|
| `data/processed/` | 4.3 GB | full-tree observables |
| `data/processed_subsampled/` | 551 MB | subsampled observables (3 108 arrays) |

## Running

```bash
julia --project=. -t auto code/1_data_generation/3_processing/<full_trees|subsampled_trees|checks>/<script>.jl
```

Files sit two levels below `1_data_generation/`, so the `ROOT` header uses four
`dirname`s.

## Resumability

- `full_trees/process_sel1.jl` and `full_trees/process_sel2.jl` are **resumable**
  per shard: a shard whose every quantity file already exists is skipped. Delete
  a processed file to force a redo.
- `full_trees/process_neutral.jl` has **no such guard**. `process_file`
  serializes unconditionally, so re-running rewrites all 433 MB of
  `data/processed/neutral/`. Do not run it to "check" something.
- Every script under `subsampled_trees/` is resumable *and* uses
  `serialize_atomic` — an interrupted run cannot leave a truncated file that the
  next run silently treats as done.

## Dependencies between scripts

**`subsampled_trees/process_sel{1,2}_subsampled.jl` require
`full_trees/process_sel{1,2}.jl` to have run first.** The subsampled processors
copy the per-simulation `injection` / `n_restarts` scalars across from
`data/processed/` rather than re-deriving them from the raw trees — re-deriving
would mean a second full pass over several GB of trees for a handful of numbers
per simulation. Both scripts hard-error, naming the required script, if those
inputs are missing or misaligned.
