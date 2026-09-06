# Inference plots

Scripts that test parameter inference against simulated data with known truth,
and later against real data.

**Empty for now, deliberately.** The inference framework lives in
`EvoTracer.jl` (separate repo, not yet implemented), and this directory holds the
study-side drivers that call it: parameter grids over selection strength and clone
birth time, sequencing-depth and sample-size sweeps, and the figures comparing
recovered parameters against the truth recorded at simulation time.

What belongs here, when it exists:

- drivers that call `EvoTracer.jl` estimators over the shards under
  `data/processed/` and `data/processed_subsampled/`
- identifiability sweeps and the figures summarising them
- comparisons of several mutation-rate estimators on the same data, including the
  Markovian estimator that is expected to fail under non-exponential timing

What does not belong here: the estimators themselves (they go in `EvoTracer.jl`,
so they can run on real data with no simulator in the dependency chain), and
anything comparing simulations to closed-form theory (that is `../theory_plots/`).

Scripts here follow the same header as every other script under `code/` — see
`code/paths.jl`.
