# `selection_1/` — exactly one driver, swept in strength and injection size

One fitness-enhancing driver mutation per simulation. Selection acts through
`f ← 1 + s`, applied to **one daughter** of the division that first takes the
population from `N_critic` to `N_critic + 1` alive cells. Uses the
`on_division` hook (MutationLoadDynamics.jl ≥ v0.2.0), which fires after the
division but *before* either daughter is scheduled — so the driver's own first
division already runs at the boosted rate, not just its descendants'.

The neutral `ν = 2.0` mutation channel is untouched (the driver does *not*
increment `mutations`), so the SFS stays directly comparable to a paired file
under `neutral/`. The driver clone is identified by `fitness > 1`; its size is
recorded in the injection record rather than in `mutations`.

Types: `X_helpers/types_selection1.jl` — `Sel1Params(b, d, gamma_shape, nu,
N_target, model, s, N_critic, rep)`, `Sel1Injection(t_inject, N_at_inject,
driver_cell_id, driver_clone_size, n_attempts, seed)`, `Sel1SimResult(trajectory,
tree_root, params, injection)`. Machinery: `X_helpers/selection1.jl`.

## Scripts and grid

Three timing-model scripts share the same three-axis grid:

- `s ∈ 0.1:0.1:2.0` — 20 driver strengths.
- `N_critic` — 10 log-spaced injection sizes from 1 to `N_target ÷ 2`
  (`ncritic_grid`). Log-spaced because `N(t) ≈ e^{rt}` in the exponential-growth
  phase, so log-spaced sizes are roughly equally-spaced injection *times*.
- 5 replicates per `(N_critic, s)` cell.

So **50 simulations per `s`** = 10 `N_critic` × 5 reps, and one output file per
`(model, N_target, d, s)` triple.

| Script | Timing | `N_target` | `d` | Shards | Sims |
|---|---|---|---|---|---|
| `growth_sel1_deterministic.jl` | `Dirac(1/(b·f))` | 1024, 16384 | 0.0 only | 40 | 2 000 |
| `growth_sel1_gamma.jl` | `Gamma(k=5, 1/(k·b·f))` | 1000, 10000 | 0.0, 0.5, 0.9 | 120 | 6 000 |
| `growth_sel1_markov.jl` | `Exponential(1/(b·f))` | 1000, 10000 | 0.0, 0.5, 0.9 | 120 | 6 000 |

**Total: 280 shards, 14 000 simulations, 6.9 GB on disk.**

## Fixed settings

- `b = 1.0`, `ν = 2.0`, driver applied via the `on_division` hook. Same
  `1e8`-scale death sentinel at `d = 0.0`.
- `restart_on_extinction = false` **plus an explicit retry loop**. A run is
  accepted only if (a) the population reached `N_target`, (b) the driver clone is
  still alive, and (c) the tree has a unique root. `n_attempts` is recorded per
  simulation; `1/mean(n_attempts)` estimates the joint survival-and-establishment
  probability.
- **Seeding hashes every field of the design point**, not just `(s, N_critic,
  rep)`. Hashing a subset made rng streams collide (`N_critic = 1` is common to
  both `N_target` values, and would have drawn the same sequence).

## Sharding, resumability, and file naming

One `.jls` file per `s` — 50 sims each. Sharding by `s` caps the worst file at
~95 MB (vs. ~1.9 GB per full parameter set) and lets a crashed sweep restart from
the first missing shard. Each script's inner loop has an `isfile(outfile) &&
continue` guard. **Delete a shard to force it to be regenerated.**

`data/raw/selection_1/{deterministic,gamma,markov}/sel1_<model>_N<N>[_d<d>][_k<k>]_s<s>.jls`

| Model | Filename example |
|---|---|
| deterministic | `sel1_deterministic_N16384_s1.5.jls` |
| gamma | `sel1_gamma_N10000_d0.5_k5.0_s1.0.jls` |
| markov | `sel1_markov_N10000_d0.9_s1.0.jls` |

## Tests and audit

- `verify_sel1.jl` — per-simulation unit suite. Writes nothing to `data/`.
- `inventory_sel1.jl` — walks every shard, checks file counts, sims per file,
  replicate/cell uniqueness, and — the part `verify_sel1.jl` cannot — that the
  `params` recorded *inside* each file agree with what its filename claims.
  Currently: `INVENTORY OK`.

## Caveats worth reading before plotting

- **`N_critic = 1` is degenerate.** Conditioned on the driver establishing, a
  driver injected at the founder's first division *is* the surviving lineage —
  at `d = 0.9, s = 0.1` the median clone is the whole population. Read that
  column separately.
- **Deterministic replicates have zero variance in the selection dynamics.**
  Timing is Dirac, so `driver_clone_size` and `n_attempts` are identical across
  replicates; only neutral mutation counts vary. Don't put error bars on those
  columns.
- **Clone size rises with `d`.** Conditioned on establishment, the driver's
  relative advantage `(b(1+s) − d)/(b − d)` is amplified near criticality — 2× at
  `d = 0` vs. 11× at `d = 0.9` for `s = 1`. Unconditional driver loss does rise
  with `d`, and shows up in `n_attempts` instead.
- **Markov needs 2–5× more attempts than gamma** at matched parameters
  (`N = 10000, d = 0.9, s = 1.0`: 14.5 vs. 3.9 attempts) — the non-Markovian
  timing effect showing up directly in the run metadata, before any SFS is
  computed.
