# `neutral/` — no selection, three timing models

Grows populations under a purely neutral process: `ν = 2.0` neutral mutations per
daughter per division, `driver_dist = Dirac(0.0)`, `fitness_update = (f, δ) -> f`.
Fitness stays exactly 1 for every cell, so `ν = 2` acts purely as the mutation rate
that feeds the SFS and the single-cell burden distribution. This is the reference
scenario against which the two selection scenarios are read.

Types: `X_helpers/types.jl` — `SimParams(b, d, gamma_shape, nu, s, N_target, N_critic)`
and `GrowthSimResult(trajectory, tree_root, params)`. `s = 0.0` for every neutral
run; `N_critic` is `NaN`.

## Scripts and grid

Three scripts, one per timing model. Each writes one `.jls` file per grid point,
holding 200 simulations.

| Script | Timing | `N_target` | `d` | Files | Sims |
|---|---|---|---|---|---|
| `growth_neutral_deterministic.jl` | `Dirac(1/b)` | 1024, 16384 | 0.0 only | 4 | 800 |
| `growth_neutral_gamma.jl` | `Gamma(k=5, 1/(k·b))` | 1000, 10000 | 0.0, 0.5, 0.9 | 6 | 1 200 |
| `growth_neutral_markov.jl` | `Exponential(1/b)` | 1000, 10000 | 0.0, 0.5, 0.9 | 6 | 1 200 |

`N_target ∈ {1024, 16384}` for deterministic (powers of two — under Dirac timing
the population doubles in lockstep, so the run overshoots to the next `2^k`);
`N_target ∈ {1000, 10000}` for gamma and markov. `k = 5` and `b = 1` are hard-coded
in the scripts, so the "grid" is a plain double loop over `(N_target, d)`.

## Fixed settings

- `b = 1.0`, `ν = 2.0`, `driver_dist = Dirac(0.0)`, `fitness_update = (f, δ) -> f`.
- **Death is exponential/gamma with mean `1/d`** and a `1e8` scale sentinel when
  `d = 0.0` so the death process never fires (identical convention across all three
  scenarios).
- `restart_on_extinction = true` — a lineage that dies before reaching `N_target`
  is restarted from a fresh founder with the same rng.
- `trajectory_dt = 0.1` in `MeasurementSpec`; only the final snapshot
  (`[AtEnd()]`) is kept — no intermediate stats are recorded.
- One `MersenneTwister(i)` seed per simulation index (thread-safe).

## Output

`data/raw/neutral/{deterministic,gamma,markov}/neutral_<model>_N<N>[_d<d>][_k<k>].jls`

| Model | Filename example |
|---|---|
| deterministic | `neutral_deterministic_N1024.jls` |
| gamma | `neutral_gamma_N10000_d0.5_k5.0.jls` |
| markov | `neutral_markov_N10000_d0.9.jls` |

`data/raw/neutral/deterministic/` also contains `N1000` and `N10000` files from
an earlier parameterization. They are orphans: downstream stages
(`process_neutral.jl`, `subsample_neutral.jl`) read only the current `N1024` and
`N16384`.

## Not resumable

None of the three scripts checks `isfile` before writing. Re-running rewrites
every output file for that model. If you only need one grid point, edit the loop
temporarily rather than trusting the sweep to skip completed shards.
