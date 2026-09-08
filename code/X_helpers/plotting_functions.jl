# Shared plotting constants and computation helpers for MutationLoadDynamics.jl.
# Include via: include(joinpath(HELPERS, "plotting_functions.jl"))
# Requires: using Statistics in the calling script.

# ── Colors ────────────────────────────────────────────────────────────────────

# Per-model data colors (mean line, uncertainty band, single-sim markers)
const COL_MARKOV        = :teal
const COL_MARKOV_BAND   = (:teal, 0.25)
const COL_MARKOV_SINGLE = (:darkgreen, 0.55)

const COL_GAMMA         = :mediumpurple
const COL_GAMMA_BAND    = (:mediumpurple, 0.25)
const COL_GAMMA_SINGLE  = (:purple4, 0.55)

const COL_DET           = :steelblue
const COL_DET_BAND      = (:steelblue, 0.25)
const COL_DET_SINGLE    = (:navy, 0.55)

# Theory reference lines — same in all model-specific plots
const COL_THEORY_MARKOV = :black
const COL_THEORY_DET    = :firebrick

# Timing distribution
const COL_DIV   = :royalblue
const COL_DEATH = :darkorange

# ── Font sizes ─────────────────────────────────────────────────────────────────
# Bumped up for publication legibility (2026-09-08): panels are also enlarged
# in the calling scripts to match, so text doesn't crowd the enlarged fonts.

const FS_HEAD   = 20
const FS_TITLE  = 18
const FS_LABEL  = 18   # axis labels read at title size
const FS_TICK   = 14
const FS_LEGEND = 15
const FS_ANNOT  = 13

# ── scMB helpers ───────────────────────────────────────────────────────────────

# Returns (hmean, hstd, single_hist) from actual per-cell mutation count vectors.
# Each element of all_mpc is a Vector{Int} of length actual_N with one mutation
# count per cell.  single_hist is from the first simulation.
function aggregate_mpc(all_mpc::Vector{Vector{Int}}, jmax::Int)
    nsims = length(all_mpc)
    mat   = zeros(Float64, jmax + 1, nsims)
    for (i, mpc) in enumerate(all_mpc)
        for m in mpc
            0 <= m <= jmax && (mat[m + 1, i] += 1.0)
        end
    end
    hmean  = vec(mean(mat; dims = 2))
    hstd   = vec(std(mat;  dims = 2))
    single = mat[:, 1]
    return hmean, hstd, single
end

# ── SFS helpers ────────────────────────────────────────────────────────────────

# Returns (smean, sstd, ssingle) from actual SFS vectors.
# Each element of all_sfs is a Vector{Int} of length actual_N where sfs[k] =
# number of mutations found in exactly k cells.  sitefrequencyspectrum already includes
# k=1 (leaf-branch / private mutations), so no correction is needed here.
function aggregate_actual_sfs(all_sfs::Vector{Vector{Int}})
    maxN = maximum(length.(all_sfs))
    mat  = zeros(Float64, maxN, length(all_sfs))
    for (i, sfs) in enumerate(all_sfs)
        mat[1:length(sfs), i] = Float64.(sfs)
    end
    smean   = vec(mean(mat; dims = 2))
    sstd    = vec(std(mat;  dims = 2))
    ssingle = Float64.(all_sfs[1])
    return smean, sstd, ssingle
end
