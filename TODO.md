# TODO

Open work, roughly in order of value. Each item says *why* it matters, not just what
to do — several look like polish but aren't.

Context: the pipeline now has five stages (see `CLAUDE.md`). Stages 1, 2 and 3 work
on complete lineage trees; stages 1b and 2b subsample. `data/` holds ~23 GB of
regenerable output, all of it gitignored.

---

## 1. Plotting for the subsampled data (stage 3b)

**Nothing under `analysis/plots/` reads `data/processed_subsampled/` yet.** 4.75 GB of
subsampled observables exist and no figure uses them, so the subsampling work is not
yet answering the question it was built for.

The arrays were deliberately laid out as a drop-in second input: they are co-indexed
with `data/processed/` (index `i` is the same simulation on both sides), and
`aggregate_actual_sfs` in `analysis/plots/plotting_functions.jl` derives `maxN` from
the data rather than assuming it, so it already accepts a length-`n` spectrum.

Suggested first figure, since it is the one the theory now makes a sharp prediction
about: the neutral SFS at full `N` versus `n = 0.1N` and `n = 0.01N`, one panel per
timing model, with the hypergeometric projection of the full spectrum overlaid on the
sampled points. Per `theory/sfs.md` the deterministic power-of-2 spikes should be
smeared below `k ≈ 4` and resolvable above it, with roughly `log2(n/4)` spikes
surviving — worth showing, because it says the fingerprint of synchronous division is
a *sample-size* question, not a population-size one.

Follow the existing convention: a panel per timing model rather than a new script per
sample size.

## 2. Update `theory/inference.md` for sampled data

`theory/sfs.md` now shows that the sample SFS is **not** Eq. 9 with `N → n`: that
substitution underestimates the singleton bin by 3.5x at 10% sampling and 7.4x at 1%.
Any inference recipe in `inference.md` that fits the closed form to sampled data
inherits that bias, most severely through `S_1`, which is also the bin with the most
counts and therefore the most fitting weight.

Two defensible fixes to document: fit the projected prediction (apply the
hypergeometric operator to the model spectrum before comparing), or restrict the fit
to intermediate `k` where the ratio is within a few percent of 1. The first is
correct; the second is cheaper and may be adequate. Worth stating which the paper
uses.

## 3. `check_sfs_projection.jl` — keep the hypergeometric oracle

`theory/sfs.md`'s projection was validated by projecting the measured full-tree SFS
and comparing against the measured subsampled SFS, agreeing to 0.0–0.12% on total
segregating sites across all three timing models. That check was run ad hoc and
discarded.

It is worth having permanently as `analysis/processing_subsampled/check_sfs_projection.jl`:
it uses **none** of the subsampling code, so it is a genuinely independent oracle; it
runs in seconds because it reads only the processed arrays; and it simultaneously
confirms that the draw is uniform without replacement, that `sfs[k]` counts mutations
in exactly `k` sampled cells, and that the induced tree retains the right mutations.
For a pipeline whose failure mode is a subtly wrong observable rather than a crash,
that is the highest-value test available.

Implementation note: compute `binomial` terms in log space (`lgamma`) — `N = 16384`
overflows otherwise.

## 4. `check_sel1_subsampled.jl` / `check_sel2_subsampled.jl`

Only the neutral scenario has a cross-check against its full-tree arrays
(`analysis/processing_subsampled/check_neutral_subsampled.jl`). The two selection
scenarios have none, which is an asymmetry with no justification behind it.

The scenario-specific quantities worth asserting:
- **sel1** — the sampled driver clone fraction (`count(>(1.0), leaf_fitness) / n`)
  against `injection.driver_clone_size / N`. A uniform sample is an unbiased estimator
  of a clone frequency, so mean absolute deviation should be near `sqrt(0.25/n)`;
  measured 0.0042–0.0052 for `n = 1000`, correlation 0.9998. This is scenario 1's
  headline observable and the reason `leaf_fitness` is stored for it at all.
- **sel2** — `cor(mut_per_cell, leaf_fitness)` per cell, which should match the
  full-tree value (measured 0.689 at `s = 0.2`, 0.943 at `s = 0.05`; the lower figure
  is the `M = 10` fitness cap compressing the top of the range). A near-zero
  correlation is the signature of the two arrays being written out of alignment.

## 5. Atomic writes in the pre-existing stages

`analysis/helpers/subsampling.jl` gained `serialize_atomic` (write to `<path>.tmp`,
then rename) because stages 1b and 2b decide "already done" from `isfile`, so an
interrupted run would otherwise leave a truncated-but-present file that every later
run silently skips. **This was not a hypothetical — it happened during development
and needed a manual audit of 50 shards to rule out corruption.**

`analysis/processing/process_sel1.jl` and `process_sel2.jl` have the identical
pattern: an `all(isfile(...))` resumability skip and a plain `serialize`. They are
exposed to exactly the same failure, over 4.3 GB of stage-2 output. Route them
through `serialize_atomic` too.

`analysis/sim_runs/**` also uses plain `serialize`, but has no `isfile` skip — an
interrupted sim run is simply re-run — so the hazard there is only that a truncated
file looks complete to a *later reader*. Lower priority, still worth doing for the
same three lines.

## 6. `process_neutral.jl` is not resumable

`process_sel1.jl` and `process_sel2.jl` skip a shard whose outputs already exist;
`process_neutral.jl` has no such check and reprocesses everything. Only 14 shards, so
it is not painful, but the inconsistency is a trap: the subsampled stages copied the
resumability idiom from `process_sel1.jl` specifically because `process_neutral.jl`
lacks it, and anyone reading the neutral script first would draw the wrong conclusion
about the convention.

## 7. A shared grid table (`analysis/helpers/grids.jl`)

The `(model, N, d, s) -> filename` grid is now written out **seven** times: three
`analysis/processing/process_*.jl`, three `analysis/subsampling/subsample_*.jl` plus
three `analysis/processing_subsampled/process_*_subsampled.jl` (which reuse the same
literals), and once more in `analysis/subsampling/inventory_subsampled.jl`.

Every copy must interpolate filenames character-for-character identically or a shard
is silently skipped with a warning rather than erroring — and `CLAUDE.md` names that
as the hazard the repetition creates. `inventory_subsampled.jl` is the compensating
control, but it is itself a further independent copy.

A single `grids.jl` yielding `(scenario, model, stem, N_target)` tuples would remove
the class of bug. It touches the pre-existing stage-2 scripts, so it wants its own
change rather than being smuggled into other work. Note the one real subtlety:
scenario 1's `s` values **must** stay `collect(0.1:0.1:2.0)` rather than a literal
list, because range arithmetic and float-literal parsing do not always agree in the
last bits and the filenames encode the value.

## 8. Small code-quality items

Deferred from review as genuinely minor. None affect results.

- `sample_seed` (`analysis/helpers/subsampling.jl:64`) has no `::UInt64` return
  annotation and relies on `hash(::Tuple)` being `UInt64` — true on 64-bit. The
  struct field is typed, so only the suite's own `isa UInt64` assertion would break
  on a 32-bit host. One token.
- The six stage scripts each carry `using AbstractTrees` and `using Random` that they
  do not use directly (`helpers/subsampling.jl` imports them itself). Harmless, and
  uniformity across the family is arguably worth more than trimming them — listed
  only so the next reader knows it was noticed and not an oversight.
- `inventory_subsampled.jl` uses top-level code with `global` declarations where
  `analysis/sim_runs/selection_1/inventory_sel1.jl` wraps the sweep in a function,
  sidestepping Julia's soft-scope rule by construction. The existing script is the
  better pattern.
- The `serialize_atomic` test asserts the end state (destination complete, no `.tmp`
  left behind) but does not simulate a mid-write interruption, so it would also pass
  against a non-atomic implementation. Testing the real thing needs a subprocess kill;
  probably not worth it, but the gap is real.

## 9. `Manifest.toml` is not portable — set up a local registry

`Manifest.toml` embeds the absolute path
`/Users/alexanderstein/Documents/GitHub/MutationLoadDynamics.jl` because the package
is added via `Pkg.develop` rather than from a registry. `Pkg.instantiate()` on any
other machine fails on that entry until it is re-`Pkg.develop`ed to wherever the
package lives there.

**This got worse on 2026-09-04.** The repo family is now three packages —
`MutationLoadDynamics.jl`, `CopyNumberEvolution.jl`, `EvoTracer.jl` — plus this study
repo and a second study repo to come (see
`docs/superpowers/specs/2026-09-04-repo-architecture-design.md`). Every one of them is
unregistered and dev-linked, so the problem triples, and a second question appears
that dev paths cannot answer at all: *which version of the sampler drew this data, and
which CNA model made these profiles?* With a dev path, `Manifest.toml` records no
version and no tree hash, so simulation output has no recoverable provenance.

Two options:

- **A private registry via `LocalRegistry.jl`.** Study repos depend on tagged
  versions; `Manifest.toml` pins a version and a tree hash instead of a local path;
  provenance becomes recordable in the output files. Roughly half an hour of setup,
  and the moment to do it is when `MutationLoadDynamics.jl` is tagged 0.3.0 for the
  sampling work. **Recommended.**
- **Keep dev paths** and add a `setup.jl` to each study repo that re-`develop`s all
  three packages from a path prefix or environment variable. Cheaper, fixes the fresh
  checkout, does nothing for provenance.

## 10. Pre-existing: `analysis/plots/old/`

Superseded by `analysis/plots/neutral/` and kept for reference. If the reference value
has expired, deleting it removes a directory that new work should not extend — the
kind of thing that is only ever obvious to the author.
