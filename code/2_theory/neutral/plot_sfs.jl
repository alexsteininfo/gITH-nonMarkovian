using Pkg
const ROOT = dirname(dirname(dirname(@__DIR__)))
isfile(joinpath(ROOT, "Project.toml")) ||
    error("ROOT = $ROOT has no Project.toml — was this script moved?")
Pkg.activate(ROOT)
include(joinpath(ROOT, "code", "paths.jl"))

using MutationLoadDynamics
using Distributions
using Random
using Serialization
using CairoMakie
using Statistics

include(joinpath(HELPERS, "types.jl"))
include(joinpath(HELPERS, "theory.jl"))
include(joinpath(HELPERS, "plotting_functions.jl"))

# ── Parameters ────────────────────────────────────────────────────────────────

const MU       = 2.0   # neutral mutations per division
const k        = 5.0
const N_VALUES = [1_000, 10_000]
const D_VALUES = [0.0, 0.5, 0.9]

const PROC = joinpath(DATA,    "processed", "neutral")
const OUT  = joinpath(FIGURES, "2_theory", "neutral")

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

# ── Figure factory ────────────────────────────────────────────────────────────

function make_sfs_figure(;
    model,
    d_values,
    col_mean, col_band, col_single,
    fig_title,
    outpath,
    n_values = N_VALUES,
)
    nrows = length(n_values)
    ncols = length(d_values)

    fig = Figure(size = (max(ncols, 2) * 380, nrows * 340 + 120), figure_padding = 18)

    Label(fig[1, 1:ncols];
        text      = fig_title,
        fontsize  = FS_HEAD, font = :bold,
        tellwidth = false, padding = (0, 0, 6, 0))

    for (ni, N) in enumerate(n_values), (di, d) in enumerate(d_values)
        fpath = proc_path(model, N, d, "sfs")
        if !isfile(fpath)
            ax = Axis(fig[ni + 1, di];
                title = ncols > 1 ? "N = $N,  d = $d" : "N = $N",
                titlesize = FS_TITLE, xscale = log10, yscale = log10, spinewidth = 1.3)
            text!(ax, 0.5, 0.5; text = "data not found",
                  align = (:center, :center), space = :relative,
                  fontsize = FS_LABEL, color = :grey60)
            println("  Missing: $fpath")
            continue
        end

        all_sfs  = deserialize(fpath)::Vector{Vector{Int}}
        actual_N = length(all_sfs[1])
        println("  $model N_target=$N N_actual=$actual_N d=$d: $(length(all_sfs)) sims")

        panel_title = ncols > 1 ? "N = $actual_N,  d = $d" : "N = $actual_N"
        ax = Axis(fig[ni + 1, di];
            title          = panel_title,
            titlesize      = FS_TITLE,
            xlabel         = ni == nrows ? "frequency k (cells)" : "",
            ylabel         = di == 1     ? "Sₖ (mutations)"      : "",
            xlabelsize     = FS_LABEL, ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,  yticklabelsize = FS_TICK,
            xscale         = log10,    yscale         = log10,
            limits         = (nothing, (0.5, nothing)),
            xgridvisible   = true,     ygridvisible   = true,
            xgridcolor     = (:grey80, 0.5), ygridcolor = (:grey80, 0.5),
            spinewidth     = 1.3)

        smean, sstd, ssingle = aggregate_actual_sfs(all_sfs)

        valid_band   = findall(smean .> 0)
        valid_single = findall(ssingle .> 0)

        if !isempty(valid_band)
            lo = max.(1e-3, smean[valid_band] .- sstd[valid_band])
            hi = smean[valid_band] .+ sstd[valid_band]
            band!(ax, Float64.(valid_band), lo, hi; color = col_band)
            lines!(ax, Float64.(valid_band), smean[valid_band];
                   color = col_mean, linewidth = 4.0)
        end

        if !isempty(valid_single)
            scatter!(ax, Float64.(valid_single), ssingle[valid_single];
                     color = col_single, markersize = 6)
        end

        # Markovian theory (Eq. 9, Stein & Werner 2025) — shown for all models.
        params  = deserialize(proc_path(model, N, d, "params"))::Vector{SimParams}
        p       = params[1]
        k_range = 1:actual_N
        theory  = [predict_sfs_theory(k, Float64(p.b), Float64(p.d), actual_N, MU)
                   for k in k_range]
        valid_th = findall(theory .> 0)
        if !isempty(valid_th)
            lines!(ax, Float64.(k_range[valid_th]), theory[valid_th];
                   color = COL_THEORY_MARKOV, linewidth = 2.5, linestyle = :dash)
        end

        # Deterministic theory — only for the perfect binary tree.
        # E[SFS[k]] = N·ν/k for k ∈ {1,2,4,…,N/2}; all other bins are zero.
        if model == "deterministic"
            L       = round(Int, log2(actual_N))
            k_det   = [1 << j for j in 0:(L - 1)]
            sfs_det = [predict_sfs_deterministic(k, actual_N, MU) for k in k_det]
            lines!(ax, Float64.(k_det), sfs_det;
                   color = COL_THEORY_DET, linewidth = 5.0, linestyle = :dot)
        end
    end

    rowgap!(fig.layout, 1, 4.0)

    leg_elems = LegendElement[
        PolyElement(color = col_band),
        LineElement(color = col_mean, linewidth = 4.0),
        MarkerElement(color = col_single, marker = :circle, markersize = 8),
        LineElement(color = COL_THEORY_MARKOV, linewidth = 2.5, linestyle = :dash),
    ]
    leg_labels = ["mean ± std", "mean", "single simulation", "Markovian theory (Eq. 9)"]
    if model == "deterministic"
        push!(leg_elems,  LineElement(color = COL_THEORY_DET, linewidth = 5.0, linestyle = :dot))
        push!(leg_labels, "Deterministic theory")
    end
    Legend(fig[nrows + 2, 1:ncols], leg_elems, leg_labels;
        orientation  = :horizontal, labelsize    = FS_LEGEND,
        framevisible = false,       tellwidth    = false,
        patchsize    = (28.0, 14.0), nbanks       = ncols > 1 ? 1 : 2)

    mkpath(dirname(outpath))
    save(outpath, fig)
    println("→ $outpath\n")
end

# ── Generate figures ──────────────────────────────────────────────────────────

println("\n── Markov SFS ──────────────────────────────────────────────────────")
make_sfs_figure(
    model      = "markov",
    d_values   = D_VALUES,
    col_mean   = COL_MARKOV,
    col_band   = COL_MARKOV_BAND,
    col_single = COL_MARKOV_SINGLE,
    fig_title  = "Site frequency spectrum — Markov model (μ = $MU)",
    outpath    = joinpath(OUT, "markov", "sfs.png"),
)

println("── Gamma SFS ───────────────────────────────────────────────────────")
make_sfs_figure(
    model      = "gamma",
    d_values   = D_VALUES,
    col_mean   = COL_GAMMA,
    col_band   = COL_GAMMA_BAND,
    col_single = COL_GAMMA_SINGLE,
    fig_title  = "Site frequency spectrum — Gamma model (k = $k, μ = $MU)",
    outpath    = joinpath(OUT, "gamma", "sfs.png"),
)

println("── Deterministic SFS ───────────────────────────────────────────────")
make_sfs_figure(
    model      = "deterministic",
    d_values   = [0.0],
    col_mean   = COL_DET,
    col_band   = COL_DET_BAND,
    col_single = COL_DET_SINGLE,
    fig_title  = "Site frequency spectrum — Deterministic model (μ = $MU)",
    outpath    = joinpath(OUT, "deterministic", "sfs.png"),
    n_values   = [1_024, 16_384],
)

println("All done.")
