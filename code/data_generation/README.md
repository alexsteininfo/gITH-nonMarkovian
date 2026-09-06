# Simulation runs

Stage 1 of the pipeline: grow populations, keep the full lineage tree, serialize to
`data/raw/`. Stage 2 is `code/data_generation/processing/`, stage 3 is `code/theory_plots/`.

Everything here runs through **MutationLoadDynamics.jl** — unregistered, under active
development, and `Pkg.develop`ed from a local path (`../MutationLoadDynamics.jl` alongside
this repo; the absolute path is baked into `Manifest.toml`). `selection_1/` requires
**v0.2.0 or later** for its `on_division` hook. If a sim script breaks after an unrelated
change, check that repo for API drift before assuming the bug is local. Run from the
**repo root**:

```bash
julia --project=. -t auto code/data_generation/<scenario>/<script>.jl
```

`-t auto` matters — every script parallelises with `Threads.@threads` and does nothing
useful single-threaded. The explicit `--project` does not rescue a wrong `Pkg.activate`
inside a script, which computes the repo root by climbing `dirname(@__DIR__)`; from a
`data_generation/<scenario>/` subdirectory that is **three** `dirname` calls to the root, **two**
to `code/`, and **three** `".."` segments to `data/`.

`data/` is gitignored and fully regenerable from these scripts.

## What has been run

| Scenario | Timing models | Files | Sims | Disk | Audit |
|---|---|---|---|---|---|
| `neutral/` | deterministic, gamma, markov | 16 | 3 200 | 1.8 GB | — |
| `selection_1/` | deterministic, gamma, markov | 280 | 14 000 | 6.9 GB | `INVENTORY OK` |
| `selection_2/` | deterministic, gamma, markov | 56 | 11 200 | 4.9 GB | `INVENTORY OK` |
| `selection_old/` | — | — | — | 733 MB | superseded |

## Shared conventions

Held identical across all three live scenarios, so files pair up across scenarios:

- `b = 1.0` — baseline birth rate at fitness `f = 1`; mean division time `1/(b·f)`.
- **Three timing models**, the recurring axis of comparison. Division waiting time for a
  cell of fitness `f`:
  - `deterministic` — `Dirac(1/(b·f))`, the `k → ∞` limit. `gamma_shape` recorded as `Inf`.
  - `gamma` — `Gamma(k, 1/(k·b·f))` with `k = 5` (CV ≈ 0.45, cancer-cell regime).
  - `markov` — `Exponential(1/(b·f))`, i.e. the `k = 1` Gamma / standard Gillespie
    birth–death model. This is the reference the other two are read against.
    `gamma_shape` recorded as `1.0`.
- **Death is fitness-independent** in every scenario. `d = 0` uses a `1e8` scale sentinel
  so death never fires; deterministic runs are `d = 0` only.
- **File naming encodes the parameters**, one file per parameter set, and omits fields that
  do not apply — so a gamma filename carries `_d` and `_k`, markov carries `_d` only, and
  deterministic carries neither. Match this for any new sweep.
- **Each scenario owns its structs** (`code/helpers/types*.jl`), `include`d rather than
  imported so independent scripts see identical definitions. They are deliberately *not*
  shared: Julia's `Serialization` resolves struct layout at deserialize time, so widening
  one struct would make every `.jls` already written against it unreadable.

## `neutral/` — no selection

Baseline. `ν = 2.0` neutral mutations per daughter per division, `driver_dist = Dirac(0.0)`
and `fitness_update = (f, δ) -> f`, so fitness stays exactly 1 and `ν` acts purely as the
mutation rate feeding the SFS and single-cell burden.

| | |
|---|---|
| Scripts | `growth_neutral_{deterministic,gamma,markov}.jl` |
| Types | `helpers/types.jl` (`SimParams`, `GrowthSimResult`) |
| Grid | gamma/markov: `N ∈ {1000, 10000}` × `d ∈ {0, 0.5, 0.9}`; deterministic: `N ∈ {1024, 16384}`, `d = 0` |
| Sims | 200 per file |
| Extinction | `restart_on_extinction = true` |

`data/raw/neutral/deterministic/` also holds `N1000` and `N10000` files from an earlier
parameterisation. `process_neutral.jl` reads only `N1024` and `N16384`; the other two are
orphans.

## `selection_1/` — exactly one driver, swept in strength and timing

One fitness-enhancing mutation per simulation, `f ← 1 + s`, applied to one daughter of the
division that first takes the population from `N_critic` to `N_critic + 1` cells. Uses the
`on_division` hook added in MutationLoadDynamics v0.2.0, which fires after the division but
*before* either daughter is scheduled — so the driver's own first division already runs at
the boosted rate, not just its descendants'.

The neutral `ν = 2.0` channel is untouched, so the driver is a separate channel and the SFS
stays directly comparable to `neutral/`. The driver does **not** increment `mutations`; the
clone is identified by `fitness > 1`.

| | |
|---|---|
| Scripts | `growth_sel1_{deterministic,gamma,markov}.jl`, `verify_sel1.jl`, `inventory_sel1.jl` |
| Types | `helpers/types_selection1.jl`, machinery in `helpers/selection1.jl` |
| Grid | `s ∈ 0.1:0.1:2.0` (20) × `N_critic` (10, log-spaced `1 → N_target/2`) × 5 reps = 1000 per parameter set |
| Sharded | one file per `s`, 50 sims each (10 `N_critic` × 5 reps) |
| Extinction | `restart_on_extinction = false` + explicit retry loop |

`N_critic` is log-spaced because `N(t) ≈ e^{rt}` during exponential growth, so log-spaced
sizes give roughly equally spaced injection *times*. A run is accepted only if the
population reached `N_target`, the driver clone is still alive, and the tree has a unique
root; `n_attempts` is recorded, and `1/E[attempts]` estimates the joint
survival-and-establishment probability.

Resumable — a shard that already exists is skipped, so delete it to force regeneration.

## `selection_2/` — every mutation is a driver, capped

No injection and no retry loop. Each daughter draws `j ~ Poisson(ν)` mutations; each draws
`X ~ Gamma(a, 1/a)` (mean 1) and updates

```
f ← min( f · (1 + s·X·(1 − f/M)),  M )
```

so `s` sets the mean per-mutation effect and `a` sets only its spread (`CV = 1/√a`). The
logistic factor gives diminishing returns; the outer `min` enforces the cap. Fitness is
non-decreasing — no deleterious mutations. At `s = 0` the rule reduces to `f ← f` exactly.

| | |
|---|---|
| Scripts | `growth_sel2_{deterministic,gamma,markov}.jl`, `verify_sel2.jl`, `inventory_sel2.jl` |
| Types | `helpers/types_selection2.jl`, machinery in `helpers/selection2.jl` |
| Grid | `s ∈ {0.05, 0.10, 0.15, 0.20}`; gamma/markov `N ∈ {1000, 10000}` × `d ∈ {0, 0.5, 0.9}`; deterministic `N ∈ {1024, 16384}` |
| Fixed | `ν = 1.0`, `M = 10.0`, `effect_shape a = 2.0` (CV ≈ 0.71), `trajectory_dt = 0.02` |
| Sims | 200 per file |
| Extinction | `restart_on_extinction = true`, `n_restarts` recorded per sim |

There is **no separate neutral channel**: the same mutations carrying the fitness effects
are the ones `sitefrequencyspectrum` and `mutations_per_cell` count. `ν = 1.0` here, not the neutral
runs' 2.0 — at `ν = 2` lineages pile onto the cap and fitness *differences* collapse.
`trajectory_dt = 0.02` rather than 0.1 because higher fitness reaches `N_target` much sooner
in simulation time.

## Auditing a completed sweep

```bash
julia --project=. code/data_generation/selection_1/inventory_sel1.jl
julia --project=. code/data_generation/selection_2/inventory_sel2.jl
```

Each walks every shard and checks file counts, sims per file, replicate/cell uniqueness, and
— the part the per-script checks cannot do — that the `params` recorded *inside* each file
agree with what its filename claims. Exit 0 and `INVENTORY OK`, or a deduplicated problem
list and exit 1. Both currently pass. `verify_sel{1,2}.jl` are the unit-level suites and
write nothing to `data/`.

## Caveats worth knowing before plotting

- **`selection_2` is near-saturated at the top two `s` values.** Mean final fitness against
  the cap `M = 10`, at `N = 10000, d = 0.5`: gamma 2.75 / 6.20 / 8.20 / 9.10 and markov
  4.50 / 7.86 / 9.06 / 9.57 for `s = 0.05 / 0.10 / 0.15 / 0.20`. The `ν = 1.0` choice was
  made expecting the worst corner near `f ≈ 6.1`; the observed values are much closer to the
  cap, so `s = 0.15` and `s = 0.20` have compressed fitness differences. `s = 0.05` and
  `s = 0.10` are cleanly graded.
- **scMB is not comparable between `neutral/` and `selection_2/`** without care: `ν = 2.0`
  versus `ν = 1.0`. SFS *shapes* still overlay after a factor-`ν` rescale; the single-cell
  burden distribution does not. A fair scMB overlay needs neutral runs at `ν = 1`.
- **`selection_1`'s `N_critic = 1` column is degenerate.** Conditioned on the driver
  establishing, a driver injected at the founder's first division *is* the surviving
  lineage — at `d = 0.9, s = 0.1` the median clone is the whole population. Read that column
  separately.
- **Deterministic replicates have zero variance in the selection dynamics.** Timing is
  Dirac, so `driver_clone_size` and `n_attempts` are identical across replicates; only
  neutral mutation counts vary. Don't put error bars on those columns.
- **Clone size rises with `d`, not falls.** Conditioned on establishment, the driver's
  relative advantage `(b(1+s) − d)/(b − d)` is amplified near criticality — 2× at `d = 0`
  versus 11× at `d = 0.9` for `s = 1`. Unconditional driver loss does rise with `d`, and
  shows up in `n_attempts` instead.
- **Markov needs 2–5× more attempts/restarts than gamma** at matched parameters
  (`selection_1`, `N = 10000, d = 0.9, s = 1.0`: 14.5 vs 3.9 attempts; `selection_2`,
  `N = 10000, d = 0.9`: 1.9–3.7 vs 1.1–1.6 restarts). Since `1/E[attempts]` is an
  establishment probability and `1/(1 + E[restarts])` a founder-survival probability, this
  is the non-Markovian timing effect showing up directly in the run metadata, before any
  SFS is computed.

## `selection_old/`

Superseded first pass at selection: four fitness models (additive, max-random, and two
multiplicative) swept over `s ∈ 0.0:0.1:0.5` at 50 sims each. Kept for reference —
`selection_2/` replaces `growth_multi_rand.jl` with a normalised effect distribution.
Don't extend it, and don't use it as a template: both its `include` and its output path are
one segment short, so it activates the wrong environment and writes to `analysis/data/`.
