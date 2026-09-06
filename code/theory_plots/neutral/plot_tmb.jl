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

const MU       = 2.0
const b        = 1.0
const k        = 5.0
const N_VALUES = [1_000, 10_000]
const D_VALUES = [0.0, 0.5, 0.9]

const PROC = joinpath(DATA,  "processed", "neutral")
const OUT  = joinpath(PLOTS, "neutral")

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

function make_tmb_figure(;
    model,
    d_values,
    col_hist, col_theory,
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
        fpath = proc_path(model, N, d, "sfs")
        if !isfile(fpath)
            ax = Axis(fig[ni + 1, di];
                title     = ncols > 1 ? "N = $N,  d = $d" : "N = $N",
                titlesize = FS_TITLE, spinewidth = 0.8)
            text!(ax, 0.5, 0.5; text = "data not found",
                  align = (:center, :center), space = :relative,
                  fontsize = 9, color = :grey60)
            println("  Missing: $fpath")
            continue
        end

        all_sfs  = deserialize(fpath)::Vector{Vector{Int}}
        actual_N = length(all_sfs[1])
        params   = deserialize(proc_path(model, N, d, "params"))::Vector{SimParams}
        p        = params[1]
        println("  $model N_target=$N N_actual=$actual_N d=$d: $(length(all_sfs)) sims")

        tmb_vals = [Float64(sum(sfs)) for sfs in all_sfs]
        theory   = predict_tmb_theory(Float64(p.b), Float64(p.d), actual_N, MU)

        panel_title = ncols > 1 ? "N = $actual_N,  d = $d" : "N = $actual_N"
        ax = Axis(fig[ni + 1, di];
            title          = panel_title,
            titlesize      = FS_TITLE,
            xlabel         = ni == nrows ? "total mutations" : "",
            ylabel         = di == 1     ? "simulations"    : "",
            xlabelsize     = FS_LABEL, ylabelsize     = FS_LABEL,
            xticklabelsize = FS_TICK,  yticklabelsize = FS_TICK,
            spinewidth     = 0.8)

        hist!(ax, tmb_vals;
              bins         = 30,
              color        = (col_hist, 0.55),
              strokewidth  = 0.5,
              strokecolor  = col_hist)

        vlines!(ax, [theory];
                color     = col_theory,
                linewidth = 2.0,
                linestyle = :dash)

        # Annotate mean and theory
        text!(ax, 0.97, 0.97;
            text  = "⟨TMB⟩ = $(round(mean(tmb_vals), digits=0))\ntheory = $(round(theory, digits=0))",
            align = (:right, :top), space = :relative,
            fontsize = FS_ANNOT, color = :grey40)
    end

    rowgap!(fig.layout, 1, 4.0)

    Legend(fig[nrows + 2, 1:ncols],
        [
            PolyElement(color = (col_hist, 0.55), strokewidth = 0.5, strokecolor = col_hist),
            LineElement(color = col_theory, linewidth = 2.0, linestyle = :dash),
        ],
        ["simulations (200)", "Markovian theory (Gunnarsson)"];
        orientation  = :horizontal, labelsize    = FS_LEGEND,
        framevisible = false,       tellwidth    = false,
        patchsize    = (16.0, 8.0))

    mkpath(dirname(outpath))
    save(outpath, fig)
    println("→ $outpath\n")
end

# ── Generate figures ──────────────────────────────────────────────────────────

println("\n── Markov TMB ──────────────────────────────────────────────────────")
make_tmb_figure(
    model      = "markov",
    d_values   = D_VALUES,
    col_hist   = COL_MARKOV,
    col_theory = COL_THEORY_MARKOV,
    fig_title  = "Total mutational burden — Markov model (μ = $MU)",
    outpath    = joinpath(OUT, "markov", "tmb.png"),
)

println("── Gamma TMB ───────────────────────────────────────────────────────")
make_tmb_figure(
    model      = "gamma",
    d_values   = D_VALUES,
    col_hist   = COL_GAMMA,
    col_theory = COL_THEORY_MARKOV,
    fig_title  = "Total mutational burden — Gamma model (k = $k, μ = $MU)",
    outpath    = joinpath(OUT, "gamma", "tmb.png"),
)

println("── Deterministic TMB ───────────────────────────────────────────────")
make_tmb_figure(
    model      = "deterministic",
    d_values   = [0.0],
    col_hist   = COL_DET,
    col_theory = COL_THEORY_MARKOV,
    fig_title  = "Total mutational burden — Deterministic model (μ = $MU)",
    outpath    = joinpath(OUT, "deterministic", "tmb.png"),
    n_values   = [1_024, 16_384],
)

println("All done.")
