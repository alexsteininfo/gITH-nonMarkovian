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

if MAX_SIMS_PER_SHARD < typemax(Int)
    println("[reduced_mode] MAX_SIMS_PER_SHARD = $MAX_SIMS_PER_SHARD (sims past this index stay `nothing`)")
end
