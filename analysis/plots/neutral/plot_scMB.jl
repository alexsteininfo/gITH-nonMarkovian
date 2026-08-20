using Pkg
Pkg.activate(dirname(dirname(dirname(@__DIR__))))

using MutationLoadDynamics
using Distributions
using Random
using Serialization
using CairoMakie
using Statistics

include(joinpath(dirname(dirname(@__DIR__)), "helpers", "types.jl"))
include(joinpath(dirname(dirname(@__DIR__)), "helpers", "theory.jl"))
include(joinpath(dirname(@__DIR__), "plotting_functions.jl"))

# ── Parameters ────────────────────────────────────────────────────────────────

const MU       = 2.0   # neutral mutations per division (embedded in simulation via ν = 2.0)
const b        = 1.0
const k        = 5.0
const N_VALUES = [1_000, 10_000]
const D_VALUES = [0.0, 0.5, 0.9]

const PROC = joinpath(@__DIR__, "..", "..", "..", "data", "processed", "neutral")
const OUT  = joinpath(@__DIR__, "..", "..", "..", "plots", "neutral")

# ── Data loader ───────────────────────────────────────────────────────────────

function proc_path(model::String, N::Int, d::Float64, quantity::String)
    if model == "deterministic"
        fname = "neutral_deterministic_N$(N).jls"
    elseif model == "markov"
        fname = "neutral_markov_N$(N)_d$(d).jls"
    else
        fname = "neutral_gamma_N$(N)_d$(d)_k$(k).jls"
    end
    return joinpath(PROC, model, quantity, fname)
end

# ── Panel draw helper ─────────────────────────────────────────────────────────

function draw_scMB_panel!(
    ax, all_mpc, _all_leaf_depths, params;
    col_mean, col_band, col_single,
)
    p         = params[1]
    actual_N  = length(all_mpc[1])
    avgen     = predict_averageGeneration(p.b, p.d, actual_N)
    avgen_det = predict_averageGeneration_deterministic(actual_N)
    jmax      = max(50, round(Int, 3 * MU * avgen))
    j_axis    = 0:jmax

    hmean, hstd, single = aggregate_mpc(all_mpc, jmax)

    # Analytic Theorem 3: Poisson(⟨L⟩) mean-field approximation for the depth distribution
    # (D_L=1 assumed). Exact in the large-N limit; small residual for d=0 (Yule, D_L=0.79).
    # See theory.md for discussion of this limitation.
    theory_m = predict_scMBdist(collect(j_axis), avgen, MU, actual_N)

    # ALTERNATIVE — empirical Theorem 3 using the actual depth distribution from simulations.
    # Avoids the D_L=1 assumption and gives exact visual agreement at all d values, at the
    # cost of depending on simulated depths rather than the fully analytic formula.
    # To activate: comment out `theory_m` above and uncomment the block below.
    # depth_pool = vcat(all_leaf_depths...)
    # freq = Dict{Int,Int}()
    # for l in depth_pool; freq[l] = get(freq, l, 0) + 1; end
    # ntot = length(depth_pool)
    # ls   = sort(collect(keys(freq)))
    # ps   = [freq[l] / ntot for l in ls]
    # theory_m = Float64[
    #     actual_N * sum(ps[i] * pdf(Poisson(MU * ls[i]), j) for i in eachindex(ls))
    #     for j in j_axis
    # ]

    theory_d = predict_scMBdist_deterministic(collect(j_axis), avgen_det, MU, actual_N)

    js = Float64.(j_axis)
    band!(ax, js, max.(0.0, hmean .- hstd), hmean .+ hstd; color = col_band)
    scatter!(ax, js, single;       color = col_single, markersize = 2.5)
    lines!(ax, js, hmean;          color = col_mean,          linewidth = 1.8)
    lines!(ax, js, theory_m;       color = COL_THEORY_MARKOV, linewidth = 1.5, linestyle = :dash)
    lines!(ax, js, theory_d;       color = COL_THEORY_DET,    linewidth = 1.5, linestyle = :dot)

    text!(ax, 0.97, 0.97;
        text  = "⟨gen⟩ ≈ $(round(avgen, digits=1))",
        align = (:right, :top), space = :relative,
        fontsize = FS_ANNOT, color = :grey50)
end

# ── Figure factory ────────────────────────────────────────────────────────────

function make_scMB_figure(;
    model,       # "markov" | "gamma" | "deterministic"
    d_values,    # e.g. [0.0, 0.5, 0.9] or [0.0]
    col_mean, col_band, col_single,
    fig_title,
    outpath,
    n_values = N_VALUES,
)
    nrows = length(n_values)
    ncols = length(d_values)

    fig = Figure(size = (max(ncols, 2) * 300, nrows * 260 + 90), figure_padding = 14)

    Label(fig[1, 1:ncols];
        text      = fig_title,
        fontsize  = FS_HEAD, font = :bold,
        tellwidth = false, padding = (0, 0, 6, 0))

    for (ni, N) in enumerate(n_values), (di, d) in enumerate(d_values)
        fpath = proc_path(model, N, d, "mut_per_cell")
        if !isfile(fpath)
            ax = Axis(fig[ni + 1, di];
                title = ncols > 1 ? "N = $N,  d = $d" : "N = $N",
                titlesize = FS_TITLE, spinewidth = 0.8)
            text!(ax, 0.5, 0.5; text = "data not found",
                  align = (:center, :center), space = :relative,
                  fontsize = 9, color = :grey60)
            println("  Missing: $fpath")
            continue
        end

        all_mpc         = deserialize(fpath)::Vector{Vector{Int}}
        all_leaf_depths = deserialize(proc_path(model, N, d, "leaf_depths"))::Vector{Vector{Int}}
        params          = deserialize(proc_path(model, N, d, "params"))::Vector{SimParams}
        actual_N        = length(all_mpc[1])
        println("  $model N_target=$N N_actual=$actual_N d=$d: $(length(all_mpc)) sims")

        panel_title = ncols > 1 ? "N = $actual_N,  d = $d" : "N = $actual_N"
        ax = Axis(fig[ni + 1, di];
            title          = panel_title,
            titlesize      = FS_TITLE,
            xlabel         = ni == nrows ? "mutations per cell, j" : "",
            ylabel         = di == 1     ? "Mⱼ (cells)"           : "",
            xlabelsize     = FS_LABEL, ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,  yticklabelsize = FS_TICK,
            spinewidth     = 0.8)

        draw_scMB_panel!(ax, all_mpc, all_leaf_depths, params;
                         col_mean, col_band, col_single)
    end

    rowgap!(fig.layout, 1, 4.0)

    Legend(fig[nrows + 2, 1:ncols],
        [
            PolyElement(color = col_band),
            LineElement(color = col_mean,          linewidth = 1.8),
            MarkerElement(color = col_single, marker = :circle, markersize = 5),
            LineElement(color = COL_THEORY_MARKOV, linewidth = 1.5, linestyle = :dash),
            LineElement(color = COL_THEORY_DET,    linewidth = 1.5, linestyle = :dot),
        ],
        ["mean ± std", "mean", "single simulation",
         "Markovian theory (Theorem 3)", "deterministic theory"];
        orientation  = :horizontal, labelsize    = FS_LEGEND,
        framevisible = false,       tellwidth    = false,
        patchsize    = (16.0, 8.0))

    mkpath(dirname(outpath))
    save(outpath, fig)
    println("→ $outpath\n")
end

# ── Generate figures ──────────────────────────────────────────────────────────

println("\n── Markov scMB ─────────────────────────────────────────────────────")
make_scMB_figure(
    model       = "markov",
    d_values    = D_VALUES,
    col_mean    = COL_MARKOV,
    col_band    = COL_MARKOV_BAND,
    col_single  = COL_MARKOV_SINGLE,
    fig_title   = "scMB — Markov model (Exponential timing, μ = $MU)",
    outpath     = joinpath(OUT, "markov", "scMB.png"),
)

println("── Gamma scMB ──────────────────────────────────────────────────────")
make_scMB_figure(
    model       = "gamma",
    d_values    = D_VALUES,
    col_mean    = COL_GAMMA,
    col_band    = COL_GAMMA_BAND,
    col_single  = COL_GAMMA_SINGLE,
    fig_title   = "scMB — Gamma model (k = $k timing, μ = $MU)",
    outpath     = joinpath(OUT, "gamma", "scMB.png"),
)

println("── Deterministic scMB ──────────────────────────────────────────────")
make_scMB_figure(
    model       = "deterministic",
    d_values    = [0.0],
    col_mean    = COL_DET,
    col_band    = COL_DET_BAND,
    col_single  = COL_DET_SINGLE,
    fig_title   = "scMB — Deterministic model (Dirac timing, μ = $MU)",
    outpath     = joinpath(OUT, "deterministic", "scMB.png"),
    n_values    = [1_024, 16_384],
)

println("All done.")
