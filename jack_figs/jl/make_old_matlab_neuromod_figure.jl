using CSV
using CairoMakie
using DataFrames
using Statistics

include(joinpath(@__DIR__, "old_matlab_neuromod_model.jl"))

const OUTPUT_DIR = joinpath(@__DIR__, "outputs")
const BIO_FIG2C_TSV = joinpath(dirname(@__DIR__), "akira", "xlsx_cache", "Fig_2C.tsv")
const OLD_MATLAB_FIGURE = joinpath(OUTPUT_DIR, "old_matlab_neuromod_figure.png")
const OLD_MATLAB_TRACES = joinpath(OUTPUT_DIR, "old_matlab_neuromod_traces.csv")
const OLD_MATLAB_POINTS = joinpath(OUTPUT_DIR, "old_matlab_neuromod_points.csv")

function load_biological_fig2c()
    df = DataFrame(CSV.File(BIO_FIG2C_TSV; delim = '\t', normalizenames = false))
    x = Float64[]
    y = Float64[]
    for row in eachrow(df)
        xv, yv = row[1], row[2]
        if ismissing(xv) || ismissing(yv) || xv == "" || yv == ""
            continue
        end
        push!(x, xv isa AbstractString ? parse(Float64, xv) : Float64(xv))
        push!(y, yv isa AbstractString ? parse(Float64, yv) : Float64(yv))
    end
    return x, y
end

function finite_xy(x, y)
    xf = Float64.(x)
    yf = Float64.(y)
    keep = isfinite.(xf) .& isfinite.(yf)
    return xf[keep], yf[keep]
end

function trend_vectors(x, y)
    x, y = finite_xy(x, y)
    length(x) < 2 && return Float64[], Float64[]
    A = hcat(ones(length(x)), x)
    b0, b1 = A \ y
    xs = collect(range(minimum(x), maximum(x); length = 100))
    return xs, b0 .+ b1 .* xs
end

function plot_old_trace!(ax, trace::DataFrame, p::OldMatlabNeuromodParams; title)
    sc1 = scale1(p)
    sci = scalei(p)
    sc3 = scale3(p)
    offsets = (si2l = 700.0, si2r = 570.0, si3l = 420.0, si3r = 270.0, si1 = 130.0)
    synscale = 30.0
    t = trace.time_s

    lines!(ax, t, trace.V1 .+ offsets.si2l; color = RGBf(0.0, 0.0, 0.55), linewidth = 1.1)
    lines!(ax, t, trace.s1 ./ sci .* synscale .+ offsets.si2l .- 80.0; color = RGBf(0.0, 0.0, 0.55), linewidth = 0.9)
    lines!(ax, t, trace.V2 .+ offsets.si2r; color = RGBf(0.0, 0.0, 0.95), linewidth = 1.1)
    lines!(ax, t, trace.s2 ./ sci .* synscale .+ offsets.si2r .- 85.0; color = RGBf(0.0, 0.0, 0.95), linewidth = 0.9)

    lines!(ax, t, trace.V3 .+ offsets.si3l; color = RGBf(0.58, 0.0, 0.0), linewidth = 1.1)
    lines!(ax, t, trace.s34 ./ sc3 .* synscale .+ offsets.si3l .- 85.0; color = RGBf(0.58, 0.0, 0.0), linewidth = 0.9)
    lines!(ax, t, trace.V4 .+ offsets.si3r; color = RGBf(0.9, 0.0, 0.0), linewidth = 1.1)
    lines!(ax, t, trace.s43 ./ sc3 .* synscale .+ offsets.si3r .- 85.0; color = RGBf(0.9, 0.0, 0.0), linewidth = 0.9)

    lines!(ax, t, trace.V0 .+ offsets.si1; color = RGBf(0.38, 0.38, 0.38), linewidth = 1.1)
    lines!(ax, t, trace.s0 ./ sc1 .* synscale .+ offsets.si1 .- 80.0; color = RGBf(0.38, 0.38, 0.38), linewidth = 0.9)

    label_x = 52.0
    text!(ax, label_x, offsets.si2l + 38; text = "Si2L", align = (:left, :center), fontsize = 13, color = RGBf(0.0, 0.0, 0.55))
    text!(ax, label_x, offsets.si2r + 38; text = "Si2R", align = (:left, :center), fontsize = 13, color = RGBf(0.0, 0.0, 0.95))
    text!(ax, label_x, offsets.si3l + 38; text = "Si3L", align = (:left, :center), fontsize = 13, color = RGBf(0.58, 0.0, 0.0))
    text!(ax, label_x, offsets.si3r + 38; text = "Si3R", align = (:left, :center), fontsize = 13, color = RGBf(0.9, 0.0, 0.0))
    text!(ax, label_x, offsets.si1 + 38; text = "Si1", align = (:left, :center), fontsize = 13, color = RGBf(0.38, 0.38, 0.38))

    x_start = p.tmax_s - 18.5
    y_start = 110.0
    lines!(ax, [x_start, x_start], [y_start, y_start + 50.0]; color = :black, linewidth = 2.0)
    lines!(ax, [x_start, x_start - 5.0], [y_start, y_start]; color = :black, linewidth = 2.0)
    text!(ax, x_start - 2.5, y_start + 15.0; text = "5 s", align = (:center, :bottom), fontsize = 12)
    text!(ax, x_start + 1.0, y_start + 25.0; text = "50 mV", rotation = pi / 2, align = (:center, :bottom), fontsize = 12)

    xlims!(ax, 50.0, p.tmax_s - 10.0)
    ylims!(ax, 40.0, 750.0)
    ax.title = title
    ax.xticksvisible = false
    ax.xticklabelsvisible = false
    ax.yticksvisible = false
    ax.yticklabelsvisible = false
    ax.xlabel = ""
    ax.ylabel = ""
    hidespines!(ax, :t, :r)
    return ax
end

function plot_comparison_panel!(ax, model_points::DataFrame)
    bio_x, bio_y = load_biological_fig2c()
    scatter!(ax, bio_x, bio_y; color = :gray45, marker = :circle, markersize = 8, strokecolor = :black, strokewidth = 0.7, label = "biology Fig. 2C")
    tx, ty = trend_vectors(bio_x, bio_y)
    lines!(ax, tx, ty; color = (:gray45, 0.8), linewidth = 2.0)

    colors = Dict("presynaptic" => RGBf(0.12, 0.42, 0.86), "postsynaptic" => RGBf(0.80, 0.22, 0.18))
    markers = Dict("presynaptic" => :circle, "postsynaptic" => :rect)
    for mode in ("presynaptic", "postsynaptic")
        sub = model_points[model_points.mode .== mode, :]
        isempty(sub) && continue
        scatter!(ax, sub.si1_frequency_hz, sub.burst_frequency_hz;
            color = colors[mode], marker = markers[mode], markersize = 14,
            strokecolor = :black, strokewidth = 1.0, label = "$(mode) trace")
    end

    ax.xlabel = "Si1 spike frequency (Hz)"
    ax.ylabel = "Si2 burst frequency (Hz)"
    ax.title = "Biological rate comparison"
    xlims!(ax, 0, max(6.5, maximum(bio_x) * 1.08))
    ylims!(ax, 0, max(0.34, maximum(bio_y) * 1.15))
    axislegend(ax; position = :lt, framevisible = false, labelsize = 11)
end

function build_old_matlab_neuromod_figure(; save_traces::Bool = false)
    mkpath(OUTPUT_DIR)

    presyn_params = old_matlab_presynaptic_params()
    postsyn_params = old_matlab_postsynaptic_params()
    @info "Running old MATLAB presynaptic-kinetics port"
    presyn_trace = simulate_old_matlab_neuromod(presyn_params)
    @info "Running old MATLAB postsynaptic-conductance port"
    postsyn_trace = simulate_old_matlab_neuromod(postsyn_params)

    if save_traces
        traces = vcat(
            hcat(DataFrame(mode = fill("presynaptic", nrow(presyn_trace))), presyn_trace),
            hcat(DataFrame(mode = fill("postsynaptic", nrow(postsyn_trace))), postsyn_trace),
        )
        CSV.write(OLD_MATLAB_TRACES, traces)
    end

    presyn_summary = trace_rate_summary(presyn_trace)
    postsyn_summary = trace_rate_summary(postsyn_trace)
    model_points = DataFrame(
        mode = ["presynaptic", "postsynaptic"],
        si1_frequency_hz = [presyn_summary.si1_frequency_hz, postsyn_summary.si1_frequency_hz],
        burst_frequency_hz = [presyn_summary.burst_frequency_hz, postsyn_summary.burst_frequency_hz],
        n_bursts = [presyn_summary.n_bursts, postsyn_summary.n_bursts],
    )
    CSV.write(OLD_MATLAB_POINTS, model_points)

    fig = Figure(size = (1500, 1050), backgroundcolor = :white, fontsize = 16)
    axA = Axis(fig[1, 1])
    axB = Axis(fig[1, 2])
    axC = Axis(fig[2, 1:2])

    plot_old_trace!(axA, presyn_trace, presyn_params; title = "Presynaptic kinetics (MATLAB params)")
    plot_comparison_panel!(axB, model_points)
    plot_old_trace!(axC, postsyn_trace, postsyn_params; title = "Postsynaptic conductance (MATLAB params)")

    text!(axA, 0.0, 1.03; text = "A", space = :relative, fontsize = 28, font = :bold, align = (:left, :bottom))
    text!(axB, 0.0, 1.03; text = "B", space = :relative, fontsize = 28, font = :bold, align = (:left, :bottom))
    text!(axC, 0.0, 1.03; text = "C", space = :relative, fontsize = 28, font = :bold, align = (:left, :bottom))

    colsize!(fig.layout, 1, Relative(0.62))
    rowsize!(fig.layout, 1, Relative(0.46))
    save(OLD_MATLAB_FIGURE, fig; px_per_unit = 2)
    return fig, model_points
end

if abspath(PROGRAM_FILE) == @__FILE__
    _, model_points = build_old_matlab_neuromod_figure()
    println("Saved figure: ", OLD_MATLAB_FIGURE)
    println("Trace CSV export is disabled by default; call build_old_matlab_neuromod_figure(save_traces = true) to write it.")
    println("Saved points: ", OLD_MATLAB_POINTS)
    println(model_points)
end
