using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Distributions
using Random
using Serialization
using CairoMakie
using Statistics

include(joinpath(dirname(@__DIR__), "helpers", "types.jl"))
include(joinpath(dirname(@__DIR__), "helpers", "theory.jl"))

# ── Parameters ────────────────────────────────────────────────────────────────

const MU       = 2.0          # neutral mutations per division (post-hoc Poisson)
const b        = 1.0
const k        = 5.0
const N_VALUES = [1_000, 10_000]
const D_VALUES = [0.0, 0.5, 0.9]

const RAW_DIR = joinpath(@__DIR__, "..", "..", "data", "raw", "growth_neutral")
const OUTDIR  = joinpath(@__DIR__, "..", "..", "plots", "scMB", "neutral")
mkpath(OUTDIR)

# ── Post-hoc neutral mutation sprinkling ─────────────────────────────────────
#
# Each cell division creates daughters that independently acquire Poisson(μ)
# neutral mutations. The founding root cell has 0 mutations (no parent division).
# We traverse top-down: each node passes inherited mutations down, and each child
# draws Poisson(μ) new mutations from its birth event (the parent's division).

function neutral_mpc(root::BinaryNode, μ::Float64, rng::AbstractRNG)
    muts  = Int[]
    stack = Tuple{BinaryNode, Int}[(root, 0)]
    while !isempty(stack)
        node, inherited = pop!(stack)
        has_left  = !isnothing(node.left)
        has_right = !isnothing(node.right)
        if !has_left && !has_right
            push!(muts, inherited)
        else
            has_left  && push!(stack, (node.left,  inherited + rand(rng, Poisson(μ))))
            has_right && push!(stack, (node.right, inherited + rand(rng, Poisson(μ))))
        end
    end
    return muts
end

function mb_histogram(mpc::Vector{Int}, jmax::Int)
    h = zeros(Float64, jmax + 1)
    for m in mpc
        0 <= m <= jmax && (h[m + 1] += 1.0)
    end
    return h
end

# ── Load and precompute all panel data ────────────────────────────────────────

struct NeutralPanelData
    j_axis      ::Vector{Int}
    hmean       ::Vector{Float64}
    hstd        ::Vector{Float64}
    single_hist ::Vector{Float64}
    theory      ::Vector{Float64}
    theory_det  ::Vector{Float64}
    avgen       ::Float64
    inferred_mu ::Vector{Float64}
end

panel_data = Matrix{Union{NeutralPanelData, Nothing}}(undef, length(N_VALUES), length(D_VALUES))

for (ni, N) in enumerate(N_VALUES), (di, d) in enumerate(D_VALUES)
    fpath = joinpath(RAW_DIR, "growth_neutral_N$(N)_d$(d)_k$(k).jls")
    if !isfile(fpath)
        panel_data[ni, di] = nothing
        println("Missing: $fpath")
        continue
    end

    sims = deserialize(fpath)
    println("N=$N, d=$d: $(length(sims)) sims loaded")

    avgen     = predict_averageGeneration(b, d, N)
    avgen_det = predict_averageGeneration_deterministic(b, d, N)
    jmax   = max(50, round(Int, 3 * MU * avgen))
    j_axis = 0:jmax

    valid = [sim for sim in sims if !isnothing(sim.tree_root)]
    isempty(valid) && (panel_data[ni, di] = nothing; continue)

    hist_matrix = zeros(Float64, jmax + 1, length(valid))
    inferred    = Float64[]
    single_hist = zeros(Float64, jmax + 1)
    got_single  = false

    for (si, sim) in enumerate(valid)
        rng = MersenneTwister(si + 7919 * ni + 97 * di)
        mpc = neutral_mpc(something(sim.tree_root), MU, rng)
        isempty(mpc) && continue

        h = mb_histogram(mpc, jmax)
        hist_matrix[:, si] = h
        if !got_single
            single_hist = copy(h)
            got_single  = true
        end

        mf::Vector{Float64} = Float64.(mpc)
        mu_mean = mean(mf)
        if mu_mean > 0
            push!(inferred, var(mf) / mu_mean - 1.0)
        end
    end

    hmean      = Vector{Float64}(vec(mean(hist_matrix; dims=2)))
    hstd       = Vector{Float64}(vec(std(hist_matrix;  dims=2)))
    theory     = Vector{Float64}(predict_scMBdist(collect(j_axis), avgen, MU, N))
    theory_det = Vector{Float64}(predict_scMBdist_deterministic(collect(j_axis), avgen_det, MU, N))

    panel_data[ni, di] = NeutralPanelData(
        collect(j_axis), hmean, hstd, single_hist, theory, theory_det, avgen, inferred
    )
    println("  avgen=$(round(avgen, digits=1)), jmax=$jmax, nsims=$(length(valid))")
end

nrows = length(N_VALUES)
ncols = length(D_VALUES)

# ──────────────────────────────────────────────────────────────────────────────
# Figure 1 — scMB_neutral.png
# ──────────────────────────────────────────────────────────────────────────────

const COL_BAND   = (:teal, 0.25)
const COL_MEAN   = :teal
const COL_SINGLE = (:darkgreen, 0.55)
const COL_THEORY = :black
const COL_DET    = :firebrick

fig1 = Figure(size=(ncols * 300, nrows * 260 + 90), figure_padding=14)

Label(fig1[1, 1:ncols];
    text      = "Single-cell mutational burden — neutral (μ = $MU)",
    fontsize  = 11, font = :bold, tellwidth = false, padding = (0, 0, 6, 0))

for (ni, N) in enumerate(N_VALUES), (di, d) in enumerate(D_VALUES)
    ax = Axis(fig1[ni + 1, di];
        title          = "N = $N,  d = $d",
        titlesize      = 10,
        xlabel         = ni == nrows ? "mutations per cell, j" : "",
        ylabel         = di == 1 ? "Mⱼ (cells)" : "",
        xlabelsize     = 9, ylabelsize     = 9,
        xticklabelsize = 8, yticklabelsize = 8,
        spinewidth     = 0.8,
    )

    pd = panel_data[ni, di]
    if isnothing(pd)
        text!(ax, 0.5, 0.5; text="data not found", align=(:center, :center),
              space=:relative, fontsize=9, color=:grey60)
        continue
    end

    js = Float64.(pd.j_axis)
    band!(ax, js, max.(0.0, pd.hmean .- pd.hstd), pd.hmean .+ pd.hstd; color=COL_BAND)
    scatter!(ax, js, pd.single_hist; color=COL_SINGLE, markersize=2.5)
    lines!(ax, js, pd.hmean;       color=COL_MEAN,   linewidth=1.8)
    lines!(ax, js, pd.theory;      color=COL_THEORY, linewidth=1.5, linestyle=:dash)
    lines!(ax, js, pd.theory_det;  color=COL_DET,    linewidth=1.5, linestyle=:dot)

    text!(ax, 0.97, 0.97;
        text  = "⟨gen⟩ ≈ $(round(pd.avgen, digits=1))",
        align = (:right, :top), space = :relative,
        fontsize = 7, color = :grey50)
end

rowgap!(fig1.layout, 1, 4.0)

Legend(fig1[nrows + 2, 1:ncols],
    [
        PolyElement(color=COL_BAND),
        LineElement(color=COL_MEAN, linewidth=1.8),
        MarkerElement(color=COL_SINGLE, marker=:circle, markersize=5),
        LineElement(color=COL_THEORY, linewidth=1.5, linestyle=:dash),
        LineElement(color=COL_DET,    linewidth=1.5, linestyle=:dot),
    ],
    ["mean ± std", "mean", "single simulation",
     "Markovian theory (Theorem 3)", "deterministic theory"];
    orientation=:horizontal, labelsize=8,
    framevisible=false, tellwidth=false, patchsize=(16.0, 8.0))

outfile1 = joinpath(OUTDIR, "scMB_neutral.png")
save(outfile1, fig1)
println("→ $outfile1")

# ──────────────────────────────────────────────────────────────────────────────
# Figure 2 — mutrate_neutral.png
# ──────────────────────────────────────────────────────────────────────────────

fig2 = Figure(size=(ncols * 260, nrows * 220 + 100), figure_padding=12)

Label(fig2[1, 1:ncols];
    text      = "Inferred neutral mutation rate (true μ = $MU)",
    fontsize  = 11, font = :bold, tellwidth = false, padding = (0, 0, 4, 0))

for (ni, N) in enumerate(N_VALUES), (di, d) in enumerate(D_VALUES)
    ax = Axis(fig2[ni + 1, di];
        title          = "N = $N,  d = $d",
        titlesize      = 10,
        xlabel         = ni == nrows ? "inferred μ" : "",
        ylabel         = di == 1 ? "simulations" : "",
        xlabelsize     = 9, ylabelsize     = 9,
        xticklabelsize = 8, yticklabelsize = 8,
        spinewidth     = 0.8,
    )

    pd = panel_data[ni, di]
    if isnothing(pd) || isempty(pd.inferred_mu)
        text!(ax, 0.5, 0.5; text="data not found", align=(:center, :center),
              space=:relative, fontsize=9, color=:grey60)
        continue
    end

    inf = pd.inferred_mu
    hist!(ax, inf; bins=30, color=(:slateblue, 0.7))
    vlines!(ax, [mean(inf)]; color=:black,    linewidth=1.5, linestyle=:dash)
    vlines!(ax, [MU];        color=:firebrick, linewidth=1.5)

    text!(ax, 0.97, 0.97;
        text  = "mean = $(round(mean(inf), digits=2))\nvar  = $(round(var(inf), digits=3))",
        align = (:right, :top), space = :relative,
        fontsize = 7, color = :black)
end

rowgap!(fig2.layout, 1, 2.0)

Legend(fig2[nrows + 2, 1:ncols],
    [
        PolyElement(color=(:slateblue, 0.7)),
        LineElement(color=:black,    linewidth=1.5, linestyle=:dash),
        LineElement(color=:firebrick, linewidth=1.5),
    ],
    ["distribution of inferred μ", "mean of estimates", "true μ = $MU"];
    orientation=:horizontal, labelsize=8,
    framevisible=false, tellwidth=false, patchsize=(16.0, 8.0))

outfile2 = joinpath(OUTDIR, "mutrate_neutral.png")
save(outfile2, fig2)
println("→ $outfile2")
