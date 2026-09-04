# gITH-nonMarkovian

How violating the Markov assumption in cell-division timing reshapes two classic
summary statistics of somatic evolution: the **site frequency spectrum (SFS)** and
the **single-cell mutational burden (scMB) distribution**.

## Background

Standard birth–death models of somatic tissue evolution assume exponentially
distributed (Markovian / Gillespie) division and death times. Real cell cycles are
not memoryless — they have a refractory minimum duration — so this project compares
three timing models for the same growth process:

| Model | Waiting time | Notes |
|---|---|---|
| **Deterministic** | fixed (Dirac) | perfect binary tree, no timing noise |
| **Gamma (k=5)** | Gamma(k, θ) | non-Markovian, CV ≈ 0.45 (cancer-cell regime) |
| **Markov** | Exponential (k=1) | the Gillespie limit |

The analytical predictions and estimators used throughout (SFS Eq. 9, scMB
Theorem 3, mutation-rate estimators) are derived in [`theory/`](theory/) and are
from **Stein & Werner (2025), *Genetics* 230(4): iyaf101**
([doi:10.1093/genetics/iyaf101](https://doi.org/10.1093/genetics/iyaf101)), plus
new work extending the mean-field (Markovian) approximation to the non-Markovian
and death-present regimes.

Simulations are run with
[**MutationLoadDynamics.jl**](https://github.com/alexanderstein/MutationLoadDynamics.jl),
a Julia package (by the same author, still under active development) that
simulates non-Markovian birth–death dynamics via a global min-heap event queue,
preserving the full lineage tree for exact SFS/coalescence analysis.

## Repository layout

```
theory/     Analytical derivations (model, SFS, scMB, inference) — read these first
analysis/
  sim_runs/               Scripts that run simulations via MutationLoadDynamics.jl → data/raw/
  processing/             Extract SFS / mutation burden / leaf depths → data/processed/
  subsampling/            Draw a uniform n-cell subsample of each tree → data/raw_subsampled/
  processing_subsampled/  Same quantities as processing/, recomputed on the subsamples
                          → data/processed_subsampled/
  plots/                  Compare empirical results against theory → plots/
  checks/                 Equivalence gate: package statistics/sampling vs. stored data
  helpers/                Shared types (types.jl) and theory formulas (theory.jl)
data/       Serialized (.jls) simulation output — raw trees and processed summaries
plots/      Generated figures (PNG)
```

## Setup

This repo has its own Julia environment (`Project.toml`/`Manifest.toml`) with
**`MutationLoadDynamics.jl` v0.3.0 or newer** added as a local **dev** dependency
(`Pkg.develop`), so it always tracks the working copy at
`/Users/alexanderstein/Documents/GitHub/MutationLoadDynamics.jl` rather than a
released version. The version floor is not optional: `sample_leaves`,
`mutations_per_cell`, `sitefrequencyspectrum`, `leaf_depths` and `leaf_fitness` all
live in the package now, and the processing and subsampling stages call them (the
package also provides `sample_trees`, `branch_spectrum` and
`filtered_mutations_per_cell`, promoted with the rest of the family but not yet
called anywhere in this repo). To instantiate on a fresh clone or another
machine:

```julia
using Pkg
Pkg.activate(".")
Pkg.develop(path="/path/to/MutationLoadDynamics.jl")
Pkg.instantiate()
```

Then run any script from the repo root, e.g.:

```sh
julia --project=. analysis/sim_runs/neutral/growth_neutral_markov.jl
julia --project=. analysis/processing/process_neutral.jl
julia --project=. analysis/plots/neutral/plot_sfs.jl
```

`data/` (regenerable simulation output, several GB) and `CLAUDE.md` are
git-ignored — see `.gitignore`.

## Status

- **Neutral** timing comparison (deterministic / gamma / markov, N ∈ {1e3, 1e4},
  d ∈ {0, 0.5, 0.9}) is fully run, processed, and plotted.
- **Selection** (4 fitness-update models: additive, max-random, multiplicative
  fixed/random effect — see [`theory/model.md`](theory/model.md)) has simulation
  scripts in place but `data/processed/selection/` and `plots/selection/` are
  still empty — not yet run in this repo.
- `analysis/plots/old/` is an earlier iteration of the plotting scripts, superseded
  by `analysis/plots/neutral/`.
