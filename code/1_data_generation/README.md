# Data generation

`code/1_data_generation/` covers stages 1, 2, 1b, 2b and 4 of the pipeline described in
`CLAUDE.md`, numbered by stage: grow populations and keep their full lineage trees
(`1_simulation_runs/neutral/`, `1_simulation_runs/selection_1/`,
`1_simulation_runs/selection_2/`), draw uniform `n`-cell subsamples of each tree
(`2_subsampling/`), extract observable arrays from both the full trees
(`3_processing/full_trees/`) and the subsampled trees (`3_processing/subsampled_trees/`),
and simulate copy-number alterations on the subsampled trees (`4_cn_evolution/`).
`3_processing/checks/` holds the gate that verifies this repo's local tree statistics
still agree with what `MutationLoadDynamics.jl` computes. Stage 3 (aggregation and
plotting against the closed-form theory) is `code/2_theory/` — an unrelated "stage 3"
name collision with `3_processing/`; the numbering here is purely about ordering within
`1_data_generation/`.

Everything here runs through **MutationLoadDynamics.jl** — unregistered, under active
development, and `Pkg.develop`ed from a local path (`../MutationLoadDynamics.jl` alongside
this repo; the absolute path is baked into `Manifest.toml`). `selection_1/` requires
**v0.2.0 or later** for its `on_division` hook. If a sim script breaks after an unrelated
change, check that repo for API drift before assuming the bug is local. Run from the
**repo root**:

```bash
julia --project=. -t auto code/1_data_generation/1_simulation_runs/<scenario>/<script>.jl
julia --project=. -t auto code/1_data_generation/2_subsampling/<script>.jl
julia --project=. -t auto code/1_data_generation/3_processing/<full_trees|subsampled_trees|checks>/<script>.jl
julia --project=. code/1_data_generation/4_cn_evolution/<script>.jl
```

`-t auto` matters for every script under `1_simulation_runs/`, `2_subsampling/` and
`3_processing/` — they parallelise with `Threads.@threads` and do nothing useful
single-threaded. The `4_cn_evolution/` scripts are single-threaded — `-t auto` is
harmless but unnecessary for them.

Every script opens with an identical six-line header, differing only in how many
`dirname`s climb back to the repo root: four for anything under `1_simulation_runs/`
or `3_processing/` (one scenario/kind subfolder below `data_generation/`), three for
`2_subsampling/` and `4_cn_evolution/` (files sit directly in them, same depth as
before the folder structures were reorganized):

```julia
using Pkg
const ROOT = dirname(dirname(dirname(dirname(@__DIR__))))  # 2_subsampling/ uses one less dirname
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))
```

`ROOT` is *asserted* to hold a `Project.toml` before `Pkg.activate` runs; `code/paths.jl`
then defines `DATA`, `FIGURES`, and `HELPERS` from that same `ROOT`. No script counts
directory levels to reach `data/` or `figures/` any more, so a script moved to the wrong
depth fails loudly at the assertion instead of silently activating the wrong environment
or writing to the wrong place — the recurring bug in the previous layout, and the reason
every file's `dirname` count had to be re-checked when this folder was restructured
(2026-09-07: flat `neutral/`, `selection_1/`, `selection_2/`, `processing/`,
`processing_subsampled/`, `subsampling/`, `checks/` became the numbered
`1_simulation_runs/`, `2_subsampling/`, `3_processing/` tree above).

`data/` is gitignored and fully regenerable from these scripts.

## What has been run

| Scenario | Timing models | Files | Sims | Disk | Audit |
|---|---|---|---|---|---|
| `1_simulation_runs/neutral/` | deterministic, gamma, markov | 16 | 3 200 | 1.8 GB | — |
| `1_simulation_runs/selection_1/` | deterministic, gamma, markov | 280 | 14 000 | 6.9 GB | `INVENTORY OK` |
| `1_simulation_runs/selection_2/` | deterministic, gamma, markov | 56 | 11 200 | 4.9 GB | `INVENTORY OK` |
| `selection_old/` (in `analysis/`) | — | — | — | 733 MB | superseded |

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
- **Each scenario owns its structs** (`code/X_helpers/types*.jl`), `include`d rather than
  imported so independent scripts see identical definitions. They are deliberately *not*
  shared: Julia's `Serialization` resolves struct layout at deserialize time, so widening
  one struct would make every `.jls` already written against it unreadable.

## `1_simulation_runs/neutral/` — no selection

Baseline. `ν = 2.0` neutral mutations per daughter per division, `driver_dist = Dirac(0.0)`
and `fitness_update = (f, δ) -> f`, so fitness stays exactly 1 and `ν` acts purely as the
mutation rate feeding the SFS and single-cell burden.

| | |
|---|---|
| Scripts | `growth_neutral_{deterministic,gamma,markov}.jl` |
| Types | `X_helpers/types.jl` (`SimParams`, `GrowthSimResult`) |
| Grid | gamma/markov: `N ∈ {1000, 10000}` × `d ∈ {0, 0.5, 0.9}`; deterministic: `N ∈ {1024, 16384}`, `d = 0` |
| Sims | 200 per file |
| Extinction | `restart_on_extinction = true` |

`data/raw/neutral/deterministic/` also holds `N1000` and `N10000` files from an earlier
parameterisation. `process_neutral.jl` reads only `N1024` and `N16384`; the other two are
orphans.

## `1_simulation_runs/selection_1/` — exactly one driver, swept in strength and timing

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
| Types | `X_helpers/types_selection1.jl`, machinery in `X_helpers/selection1.jl` |
| Grid | `s ∈ 0.1:0.1:2.0` (20) × `N_critic` (10, log-spaced `1 → N_target/2`) × 5 reps = 1000 per parameter set |
| Sharded | one file per `s`, 50 sims each (10 `N_critic` × 5 reps) |
| Extinction | `restart_on_extinction = false` + explicit retry loop |

`N_critic` is log-spaced because `N(t) ≈ e^{rt}` during exponential growth, so log-spaced
sizes give roughly equally spaced injection *times*. A run is accepted only if the
population reached `N_target`, the driver clone is still alive, and the tree has a unique
root; `n_attempts` is recorded, and `1/E[attempts]` estimates the joint
survival-and-establishment probability.

Resumable — a shard that already exists is skipped, so delete it to force regeneration.

## `1_simulation_runs/selection_2/` — every mutation is a driver, capped

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
| Types | `X_helpers/types_selection2.jl`, machinery in `X_helpers/selection2.jl` |
| Grid | `s ∈ {0.05, 0.10, 0.15, 0.20}`; gamma/markov `N ∈ {1000, 10000}` × `d ∈ {0, 0.5, 0.9}`; deterministic `N ∈ {1024, 16384}` |
| Fixed | `ν = 1.0`, `M = 10.0`, `effect_shape a = 2.0` (CV ≈ 0.71), `trajectory_dt = 0.02` |
| Sims | 200 per file |
| Extinction | `restart_on_extinction = true`, `n_restarts` recorded per sim |

There is **no separate neutral channel**: the same mutations carrying the fitness effects
are the ones `sitefrequencyspectrum` and `mutations_per_cell` count. `ν = 1.0` here, not the neutral
runs' 2.0 — at `ν = 2` lineages pile onto the cap and fitness *differences* collapse.
`trajectory_dt = 0.02` rather than 0.1 because higher fitness reaches `N_target` much sooner
in simulation time.

## `3_processing/full_trees/` — stage 2: full-tree observables

One script per scenario. Each deserializes the stage-1 raw trees for that scenario and
walks them to extract the arrays everything downstream consumes, one quantity per
subdirectory under `data/processed/<scenario>/`:

| | |
|---|---|
| Scripts | `process_neutral.jl`, `process_sel1.jl`, `process_sel2.jl` |
| Quantities (all three) | `params`, `mut_per_cell`, `sfs`, `leaf_depths` |
| Extra (`selection_1`) | `injection` — `Sel1Injection`: `t_inject`, `N_at_inject`, `driver_cell_id`, `driver_clone_size`, `n_attempts`, `seed` |
| Extra (`selection_2`) | `leaf_fitness`, `n_restarts` |

`process_sel1.jl` and `process_sel2.jl` skip a raw file once every one of its quantity
files already exists under `data/processed/`, so an interrupted run resumes and a forced
redo is a file deletion. **`process_neutral.jl` has no such guard** — `process_file` there
serializes unconditionally, with no per-shard `isfile` check at all — so re-running it
rewrites all 433 MB of `data/processed/neutral/`. **Never run `process_neutral.jl`, or its
port, to "check" something**; there is nothing it skips.

## `2_subsampling/` — stage 1b: uniform subsamples

For each full tree, draws `n` of its `N` leaf cells uniformly without replacement and
serializes the *induced* tree — the sampled leaves plus every ancestor of a sampled leaf,
with unary nodes retained rather than collapsed — to `data/raw_subsampled/<scenario>/`.
Retaining every ancestor is what makes a sampled cell's `mut_per_cell`/`leaf_depths` in
stage 2b identical to its full-tree value; see `theory/sfs.md` and `X_helpers/subsampling.jl`
for why, and `sample_sizes` in that file for the fixed `N_target -> [n, ...]` table (adding
a population size means adding a row there, not computing one).

| | |
|---|---|
| Scripts | `subsample_neutral.jl`, `subsample_sel1.jl`, `subsample_sel2.jl` |
| Tests | `verify_subsampling.jl` — unit suite; its strongest check draws `n = N_full` and asserts the result reproduces `data/processed/` exactly |
| Audit | `inventory_subsampled.jl` — rebuilds every expected filename from the parameter grids and compares against disk (`--deep` additionally reads every file back) |

Writes are atomic (write to `<path>.tmp`, then rename) so an interrupted run cannot leave
a truncated file that a later run mistakes for "already done".

## `3_processing/subsampled_trees/` — stage 2b: observables from the subsamples

The same quantities `full_trees/` extracts, recomputed on the subsampled trees, into
`data/processed_subsampled/<scenario>/`. `sfs` is the quantity sampling genuinely
distorts (see `theory/sfs.md`'s hypergeometric projection); `mut_per_cell` and
`leaf_depths` are unchanged from the full-tree values for the same cells, by
construction of stage 1b.

| | |
|---|---|
| Scripts | `process_neutral_subsampled.jl`, `process_sel1_subsampled.jl`, `process_sel2_subsampled.jl` |
| Cross-check | `check_neutral_subsampled.jl` — asserts the subsampled arrays against the full-tree arrays they are co-indexed with; there is no `sel1`/`sel2` equivalent yet (see `TODO.md`) |

**`process_sel1_subsampled.jl` and `process_sel2_subsampled.jl` require
`full_trees/process_sel1.jl` and `full_trees/process_sel2.jl`, respectively, to have run
first.** Each copies its scenario's per-simulation scalars (`injection` for `sel1`,
`n_restarts` for `sel2`) across from `data/processed/` rather than re-deriving them from
the raw trees — re-deriving would mean a second full pass over several GB of trees for a
handful of numbers per simulation. Both hard-error, naming the required script, if those
inputs are missing or misaligned.

## `3_processing/checks/` — cross-package equivalence gate

`check_package_equivalence.jl` is the gate that verifies this repo's own tree-statistics
helpers and sampler still agree with what `MutationLoadDynamics.jl` produces, so far
enough into the package's development that the local implementations they replaced could
be retired. It runs two testsets against real data: does the package's tree-statistics
code return exactly what `helpers/tree_analysis.jl` returns, element for element and in
the same order (`data/processed/`), and does the package's sampler reproduce the draws
already serialized under `data/raw_subsampled/`. It exits non-zero on any mismatch; a
mismatch means the package integration is wrong, never that the data should be
regenerated.

## `4_cn_evolution/` — stage 4: copy-number evolution

Simulates allele-specific copy-number alterations on the stage-1b subsampled
trees with **CopyNumberEvolution.jl** (a separate, unregistered, actively
developed package, `Pkg.develop`ed from `../CopyNumberEvolution.jl`), and writes
MEDICC2-ready output plus ground truth to `data/CN_subsampled/`. Design rationale
and every scoping decision below: `docs/superpowers/specs/2026-09-07-cn-evolution-design.md`.

Model: `CNAModel(rate = FromEdgeMutations())` on `hg38(:female)` with 1Mb bins —
every other component (target, extent, gain/loss, WGD, viability, root state) is
the package default, which already means focal-only events and no WGD. The rate
rule is exact identity: a tree edge's already-recorded mutation count becomes its
CNA count, one for one.

| | |
|---|---|
| Scripts | `cn_neutral.jl` (all 3 models), `cn_sel1.jl`, `cn_sel2.jl` (gamma only) |
| Shared driver | `code/X_helpers/cn_evolution.jl` — model constants, seed fn, `cn_evolve_sim`/`cn_evolve_file` |
| Source | `data/raw_subsampled/`, smallest sample size per shard family only (`n = 100/102/164`, never `1000`/`1638`) |
| Sims per shard | 1 (`SIMS_PER_SHARD` in `cn_evolution.jl`) — raise later; resumable, so raising it only adds work |
| Shards in scope | 158: 14 neutral (all models) + 120 sel1-gamma + 24 sel2-gamma |
| Output | `data/CN_subsampled/<scenario>/<model>/<stem>/sim<sim_index>_{cells.tsv,truth_profiles.tsv,truth_events.tsv,tree.nwk}` |
| Disk | 2.1 GB measured (1.6 GB of it `cells.tsv`) |

`selection_1`/`selection_2`'s deterministic and markov models are not yet run —
deferred for disk/runtime reasons during design, not a modelling choice. Adding
them is a `markov`/`deterministic` loop in `cn_sel1.jl`/`cn_sel2.jl`, the same
shape `subsample_sel1.jl`/`subsample_sel2.jl` already show.

No `verify_*.jl`/`inventory_*.jl` yet — premature before `SIMS_PER_SHARD` and the
1Mb resolution are validated against real usage. `medicc2` itself was not run
against the output on the machine this stage was built on (not installed); the
format was checked structurally instead (header line, autosome-only chrom column,
copy numbers within MEDICC2's 0–8 range).

## Auditing a completed sweep

```bash
julia --project=. code/1_data_generation/1_simulation_runs/selection_1/inventory_sel1.jl
julia --project=. code/1_data_generation/1_simulation_runs/selection_2/inventory_sel2.jl
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

Superseded first pass at selection (four fitness models — additive, max-random, and two
multiplicative — swept over `s ∈ 0.0:0.1:0.5` at 50 sims each). Not ported to `code/`: it
remains at `analysis/sim_runs/selection_old/` and will be deleted along with the rest of
`analysis/`. Don't use it as a template: both its `include` and its output path are one
segment short, so it activates the wrong environment and writes to `analysis/data/`.

## A note on `analysis/` references

Comments and docstrings under `code/` still name `analysis/` paths in a few dozen places
— provenance notes such as "mirrors `analysis/processing/process_sel1.jl`" and the like.
They are accurate while both trees exist, since the scripts they annotate really were
ported from those originals, and are left alone deliberately. They should be swept in the
same change that deletes `analysis/` and updates `README.md`/`CLAUDE.md`.
