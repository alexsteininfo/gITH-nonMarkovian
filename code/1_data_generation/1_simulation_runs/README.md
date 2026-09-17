# `1_simulation_runs/` — stage 1: grow populations, keep their lineage trees

Every script here calls **MutationLoadDynamics.jl** to grow a single population from
one founder cell up to a target size `N_target`, retains the entire binary lineage
tree (`BinaryNode{NonMarkovCell}`), and serializes a `Vector{<Scenario>SimResult}`
under `data/raw/<scenario>/<model>/`. Nothing here reads previous stages; this is
the source of the whole pipeline.

Three scenarios, three timing models each, run in a fully-crossed grid on
`(N_target, d, s, k)` — subfolders below.

## What is fixed across every scenario

- `b = 1.0` — baseline birth rate at fitness `f = 1`; mean division time `1/(b·f)`.
- **Three timing models**, the recurring axis of comparison. Division waiting time
  for a cell of fitness `f`:
  - `deterministic` — `Dirac(1/(b·f))`, the `k → ∞` limit. `gamma_shape = Inf`.
  - `gamma` — `Gamma(k, 1/(k·b·f))` with `k = 5` (CV ≈ 0.45, cancer-cell regime).
    `gamma_shape = 5.0`.
  - `markov` — `Exponential(1/(b·f))`, i.e. the `k = 1` Gamma / standard Gillespie
    birth–death process. `gamma_shape = 1.0`.
- **Death is fitness-independent** everywhere. `d = 0` uses a `1e8`-scale sentinel
  so death never fires; deterministic runs are `d = 0` only.
- **Filenames encode every parameter** — a gamma file carries `_d` and `_k`, markov
  carries `_d` only, deterministic carries neither `_d` nor `_k`. `neutral` files
  omit `_s`; both selection scenarios include it. Match this if you add a sweep.
- **Each scenario owns its structs** (`code/X_helpers/types*.jl`), `include`d
  rather than imported so independent scripts see identical struct layouts.
  Renaming a field of `SimParams` / `Sel1Params` / `Sel2Params` would make every
  already-written `.jls` (14 GB of raw trees alone) unreadable.
- Every script parallelizes with `Threads.@threads`, so run with `-t auto`.

## What has been run

| Scenario | Models | Files | Sims/file | Total sims | Disk |
|---|---|---|---|---|---|
| `neutral/` | deterministic, gamma, markov | 16 | 200 | 3 200 | 1.8 GB |
| `selection_1/` | deterministic, gamma, markov | 280 | 50 | 14 000 | 6.9 GB |
| `selection_2/` | deterministic, gamma, markov | 56 | 200 | 11 200 | 4.9 GB |

`inventory_sel1.jl` / `inventory_sel2.jl` audit the two selection sweeps by
walking every shard and cross-checking file counts, sims per file, replicate
uniqueness and the `params` recorded *inside* each file against what the filename
claims. Both currently pass (`INVENTORY OK`). `verify_sel{1,2}.jl` are the
per-simulation unit suites; neither writes to `data/`.

## Running

From the repo root, every script here activates the repo's environment through the
`ROOT = dirname(dirname(dirname(dirname(@__DIR__))))` header (four levels — files
sit two levels under `1_data_generation/`):

```bash
julia --project=. -t auto code/1_data_generation/1_simulation_runs/<scenario>/<script>.jl
```

Resumability differs by scenario — see each subfolder's README. In brief:

- `neutral/` — **not resumable**, `serialize` is unconditional. Re-running rewrites
  every output file for its scenario.
- `selection_1/` — resumable per shard (one file per `s`). Delete a shard to redo.
- `selection_2/` — **not resumable** at the shard level (`serialize` is
  unconditional inside the `(N, d, s)` triple loop). Resumability lives one stage
  later, in `3_processing/full_trees/process_sel2.jl`.
