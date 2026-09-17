# Env-var control for a reduced-scope run of the stage-3 drivers.
#
# `INFERENCE_MAX_SIMS` caps the number of sims processed per shard. Every driver
# uses this to bound its inner Threads.@threads loop. Sims past the cap are left
# as `nothing` in the output vector, which the figure code already skips.
#
# Default: no cap. Set `INFERENCE_MAX_SIMS=20` to smoke-test the whole pipeline
# in minutes rather than days.
#
# **Warning.** Because the output path is unchanged, running the reduced pipeline
# leaves partial `.jls` files at `data/inference/…/<stem>.jls` that the full
# pipeline's resumability check will treat as complete. `rm -rf data/inference/`
# before switching from reduced to full.

const MAX_SIMS_PER_SHARD = let s = get(ENV, "INFERENCE_MAX_SIMS", "")
    isempty(s) ? typemax(Int) : parse(Int, s)
end

# Shard-skip regexes: separate ones for full-tree and subsampled-tree drivers,
# because the same filename fragment (e.g., `_N10000_`) can mean "skip" in a
# full-tree file but "keep" in a subsampled 1%-of-N=10000 file.
#
# - `INFERENCE_SHARD_SKIP_FULL` — matched against full-tree filenames
# - `INFERENCE_SHARD_SKIP_SUB`  — matched against subsampled-tree filenames
# - `INFERENCE_SHARD_SKIP`      — legacy default for both, if the specific ones
#                                 are unset
#
# Example: keep only small-N full trees + 1% subsamples of the big-N runs:
#
#     INFERENCE_SHARD_SKIP_FULL='N10000|N16384' \
#     INFERENCE_SHARD_SKIP_SUB='_N1000_|_N1024|_n1000\.jls$|_n1638\.jls$' \
#     bash code/3_inference/run_pipeline.sh
const _SHARD_SKIP_LEGACY = get(ENV, "INFERENCE_SHARD_SKIP", "")

const SHARD_SKIP_FULL = let s = get(ENV, "INFERENCE_SHARD_SKIP_FULL", _SHARD_SKIP_LEGACY)
    isempty(s) ? nothing : Regex(s)
end

const SHARD_SKIP_SUB = let s = get(ENV, "INFERENCE_SHARD_SKIP_SUB", _SHARD_SKIP_LEGACY)
    isempty(s) ? nothing : Regex(s)
end

# Single-arg form defaults to the FULL regex — full-tree drivers keep their
# existing call site unchanged. Subsampled drivers pass `true`.
shard_should_skip(filename::AbstractString) = shard_should_skip(filename, false)
function shard_should_skip(filename::AbstractString, is_subsampled::Bool)
    r = is_subsampled ? SHARD_SKIP_SUB : SHARD_SKIP_FULL
    r !== nothing && occursin(r, filename)
end

if MAX_SIMS_PER_SHARD < typemax(Int) || SHARD_SKIP_FULL !== nothing || SHARD_SKIP_SUB !== nothing
    println("[reduced_mode] MAX_SIMS_PER_SHARD=$MAX_SIMS_PER_SHARD  SKIP_FULL=$(SHARD_SKIP_FULL)  SKIP_SUB=$(SHARD_SKIP_SUB)")
end
