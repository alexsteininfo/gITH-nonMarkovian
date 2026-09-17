# `3_processing/checks/` — cross-package equivalence gate

`check_package_equivalence.jl` is the gate that verifies this repo's stored
observables and stored subsampling draws still agree with what
**MutationLoadDynamics.jl** produces today. It exists because every tree
statistic and the leaf-sampling routine used to live locally in `analysis/`
helpers; the migration to the package (v0.3.0+) is only safe as long as this
gate keeps passing.

## What it asserts

Two testsets, both run against real data on disk:

1. **Tree statistics agree.** For every raw shard under `data/raw/`, redeserialize
   the trees, recompute `mutations_per_cell`, `sitefrequencyspectrum` and
   `leaf_depths` (and `leaf_fitness` for `selection_2`) with the package, and
   assert element-for-element equality against the arrays already stored under
   `data/processed/`. Same order, same values.
2. **Subsample draws reproduce.** For every subsampled shard under
   `data/raw_subsampled/`, re-run `sample_leaves` from the stored `seed` and
   assert the resulting leaf-id set matches the stored `sampled_ids`.

Exit 0 on complete agreement, non-zero on any mismatch. A mismatch means the
package integration is wrong; it never means the data on disk should be
regenerated.

## Running

```bash
julia --project=. -t auto code/1_data_generation/3_processing/checks/check_package_equivalence.jl
```

Reasonable to run whenever `MutationLoadDynamics.jl` is bumped locally, and
before any change to the processing scripts here.
