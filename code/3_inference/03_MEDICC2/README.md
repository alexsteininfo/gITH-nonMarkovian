# MEDICC2 sweep

Runs [MEDICC2](https://bitbucket.org/schwarzlab/medicc2) on every
`sim*_cells.tsv` under `data/CN_subsampled/`, producing a reconstructed tree,
ancestral (internal-node) copy-number profiles and copy-number events for each
simulation.

## Layout

Input (already prepared upstream):

    data/CN_subsampled/<scenario>/<model>/<param_dir>/sim<K>_cells.tsv

Output (mirrors the input, one directory per simulation):

    data/MEDICC2/treeinference/<scenario>/<model>/<param_dir>/sim<K>/sim<K>_*
    data/MEDICC2/treeinference/logs/                    # per-sim stdout+stderr

`<scenario>` is `neutral`, `selection_1` or `selection_2`; `<model>` is
`deterministic`, `gamma` or `markov`. Currently there is one simulation per
`<param_dir>` (`sim1`), 158 in total.

## Running

Intended to be run in a dedicated tmux session on this shared HPC. MEDICC2 uses
32 cores per invocation; simulations are processed one after another (do not
launch multiple copies of the runner in parallel).

    chmod +x code/3_inference/03_MEDICC2/treeinference/run_medicc2_all.sh
    ./code/3_inference/03_MEDICC2/treeinference/run_medicc2_all.sh

MEDICC2 v1.0.2 from `/srv/home/aste0033/micromamba/envs/python/bin/medicc2` is
called directly -- no `micromamba run` wrapper needed.

## Runtime

Rough guess: on the order of a few hours end-to-end (~1-3 min per simulation
×158 simulations). Not measured -- treat this as a very loose upper bound
until the first few runs report their own timings.

## Resumability

Each simulation is skipped when its `<prefix>_final_tree.new` already exists,
so an interrupted run is continued by simply re-running the script. To force a
rerun of one simulation, delete its output subdirectory. Never wipe
`data/MEDICC2/treeinference/` to "start clean" without a good reason -- the sweep is
expensive to redo.

## MEDICC2 arguments

    --input-type t          # TSV input (already the CN_subsampled format)
    --prefix sim<K>         # matches the input filename stem
    --n-cores 32
    --events                # copy-number change events, per the request
    --plot none             # no built-in plots; we produce our own downstream

No `--filter-segment-length`: the inputs are clean 1 Mb bins from simulation,
not noisy real-data CN calls, so segment-length filtering (the template's
default for the mpnst-scOmics DLP pipeline) is not needed here. Ancestral
copy-number profiles for internal tree nodes appear in
`<prefix>_final_cn_profiles.tsv` by default.
