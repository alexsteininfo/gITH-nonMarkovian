# `selection_2/` — every mutation is a driver, capped at `M`

Fully distributed selection: no injection, no retry loop. At every division each
daughter draws `j ~ Poisson(ν)` mutations, and each mutation draws
`X ~ Gamma(a, 1/a)` (mean 1) and updates

```
f ← min( f · (1 + s · X · (1 − f/M)),  M )
```

so `s` sets the *mean* per-mutation effect, `a = EFFECT_SHAPE` sets its spread
(`CV = 1/√a ≈ 0.71`), the logistic factor `(1 − f/M)` gives diminishing returns,
and the outer `min` enforces the cap. **Fitness is non-decreasing** — no
deleterious mutations. At `s = 0` the rule reduces to `f ← f` exactly.

There is **no separate neutral channel**: the same mutations that carry the
fitness effects are the ones `sitefrequencyspectrum` and `mutations_per_cell`
count. `mut_per_cell` is therefore also the driver count per cell, co-indexed with
`leaf_fitness` (both iterate `Leaves` in the same order).

Types: `X_helpers/types_selection2.jl` — `Sel2Params(b, d, gamma_shape, nu,
N_target, model, s, M, effect_shape, rep, seed)`, `Sel2SimResult(trajectory,
tree_root, params, n_restarts)`. Machinery: `X_helpers/selection2.jl`.

## Scripts and grid

Three timing-model scripts, one file per `(model, N_target, d, s)`, 200 sims each.

- `s ∈ {0.05, 0.10, 0.15, 0.20}`.
- `N_target ∈ {1000, 10000}` (gamma/markov), `{1024, 16384}` (deterministic).
- `d ∈ {0.0, 0.5, 0.9}` (gamma/markov), `0.0` only (deterministic).

| Script | Timing | `N_target` | `d` | Shards | Sims |
|---|---|---|---|---|---|
| `growth_sel2_deterministic.jl` | `Dirac(1/(b·f))` | 1024, 16384 | 0.0 only | 8 | 1 600 |
| `growth_sel2_gamma.jl` | `Gamma(k=5, 1/(k·b·f))` | 1000, 10000 | 0.0, 0.5, 0.9 | 24 | 4 800 |
| `growth_sel2_markov.jl` | `Exponential(1/(b·f))` | 1000, 10000 | 0.0, 0.5, 0.9 | 24 | 4 800 |

**Total: 56 shards, 11 200 simulations, 4.9 GB on disk.**

## Fixed settings

- `b = 1.0`, `k = 5.0` (gamma only).
- `ν = 1.0` — **not** the neutral runs' `2.0`. At `ν = 2` lineages accumulate
  ~27 mutations by `N = 10000` and the population piles onto the cap, collapsing
  fitness *differences* into an effectively neutral process with a faster clock.
  At `ν = 1` the worst corner (`s = 0.2, N = 10000`, ~14 mutations per lineage)
  reaches `f ≈ 6.1` of `M = 10` with a 5–95 percentile spread of `[4.2, 8.1]`, so
  all four `s` values stay graded (see the caveat below for observed corner-case
  values).
- `M = 10.0` — fitness cap.
- `EFFECT_SHAPE (a) = 2.0` — `CV = 1/√2 ≈ 0.71`.
- `trajectory_dt = 0.02` — not the neutral runs' `0.1`, because higher fitness
  reaches `N_target` much sooner in simulation time.
- `restart_on_extinction = true`, and `n_restarts` (extinctions before the
  attempt that reached `N_target`) is recorded per simulation.
- Seeding is per-`(model, N, d, s, rep)` via `sel2_seed`.

## Sharding, resumability, and file naming

**Not resumable at the shard level.** The inner triple-loop `serialize`s
unconditionally. Resumability lives one stage down, in
`3_processing/full_trees/process_sel2.jl`. If you need to redo one shard, either
delete-and-run or edit the loop; a full re-run rewrites all 4.9 GB.

`data/raw/selection_2/{deterministic,gamma,markov}/sel2_<model>_N<N>[_d<d>][_k<k>]_s<s>_M<M>.jls`

| Model | Filename example |
|---|---|
| deterministic | `sel2_deterministic_N16384_s0.15_M10.0.jls` |
| gamma | `sel2_gamma_N10000_d0.5_k5.0_s0.15_M10.0.jls` |
| markov | `sel2_markov_N10000_d0.9_s0.15_M10.0.jls` |

Note that `0.10` and `0.20` interpolate as `s0.1` and `s0.2` in filenames.

## Tests and audit

- `verify_sel2.jl` — per-simulation unit suite. Writes nothing to `data/`.
- `inventory_sel2.jl` — walks every shard, cross-checks file counts, sims per
  file, replicate uniqueness and in-file `params` against filename claims.
  Currently: `INVENTORY OK`.

## Caveats worth reading before plotting

- **Near-saturation at the top two `s`.** Mean final fitness against the cap
  `M = 10`, at `N = 10000, d = 0.5`:
  - gamma: `s = 0.05 → 2.75`, `0.10 → 6.20`, `0.15 → 8.20`, `0.20 → 9.10`
  - markov: `s = 0.05 → 4.50`, `0.10 → 7.86`, `0.15 → 9.06`, `0.20 → 9.57`

  The `ν = 1.0` choice was made expecting the worst corner near `f ≈ 6.1`; the
  observed values are much closer to the cap, so `s = 0.15` and `s = 0.20` have
  compressed fitness differences. `s = 0.05` and `s = 0.10` are cleanly graded.
- **scMB is not comparable to `neutral/` without care.** `ν = 1.0` here vs.
  `ν = 2.0` for `neutral/`. SFS *shapes* still overlay after a factor-`ν`
  rescale; the single-cell burden distribution does not. A fair scMB overlay
  needs neutral runs at `ν = 1`.
- **Markov needs 2–3× more restarts than gamma** at matched parameters
  (`N = 10000, d = 0.9`: 1.9–3.7 vs. 1.1–1.6). `1/(1 + mean(n_restarts))` is a
  founder-survival probability — the non-Markovian timing effect showing up in
  the metadata before any observable is computed.
