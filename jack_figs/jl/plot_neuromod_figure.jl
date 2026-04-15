using CSV
using CairoMakie
using DataFrames
using FileIO
using GLMakie

const OUTPUT_DIR = joinpath(@__DIR__, "outputs")
const PAPER_FIG5_JPG = joinpath(dirname(dirname(@__DIR__)), "paper", "fig5.jpg")
const LEGACY_CIRCUIT_SOURCE = joinpath(dirname(@__DIR__), "s1_Fig2A_drives_cpg_als_circuit.jpg")
const LEGACY_CIRCUIT_CROP = joinpath(OUTPUT_DIR, "legacy_circuit_crop2.jpg")
const FIGURE_PNG = joinpath(OUTPUT_DIR, "neuromod_master_figure.png")

const AKIRA_DIR = joinpath(dirname(@__DIR__), "akira")
const AKIRA_EXTRACTED = joinpath(AKIRA_DIR, "extracted")
const BIO_SI1_TXT = joinpath(AKIRA_EXTRACTED, "130618-03 Si1.txt")
const BIO_SI2_TXT = joinpath(AKIRA_EXTRACTED, "130618-03 Si2.txt")
const BIO_FIG2C_TSV = joinpath(AKIRA_DIR, "xlsx_cache", "Fig_2C.tsv")

const PRESYNAPTIC_POINTS_CSV = joinpath(OUTPUT_DIR, "presynaptic_scan_points.csv")
const PRESYNAPTIC_SUMMARY_CSV = joinpath(OUTPUT_DIR, "presynaptic_scan_summary.csv")
const POSTSYNAPTIC_POINTS_CSV = joinpath(OUTPUT_DIR, "postsynaptic_scan_points.csv")
const POSTSYNAPTIC_SUMMARY_CSV = joinpath(OUTPUT_DIR, "postsynaptic_scan_summary.csv")
const BIO_DRIVEN_TRACES_CSV = joinpath(OUTPUT_DIR, "bio_driven_traces.csv")
const BIO_DRIVEN_POINTS_CSV = joinpath(OUTPUT_DIR, "bio_driven_points.csv")

const MODE_COLORS = Dict(
    "presynaptic" => RGBf(0.12, 0.42, 0.86),
    "postsynaptic" => RGBf(0.80, 0.22, 0.18),
)
const MODE_MARKERS = Dict(
    "presynaptic" => :circle,
    "postsynaptic" => :rect,
)

function clean_name(x)
    replace(String(x), "_" => " ")
end

function load_two_column_txt(path::AbstractString)
    df = DataFrame(CSV.File(path; delim = '\t', normalizenames = false))
    return Float64.(df[!, 1]), Float64.(df[!, 2])
end

function load_biological_fig2c()
    df = DataFrame(CSV.File(BIO_FIG2C_TSV; delim = '\t', normalizenames = false))
    x = Float64[]
    y = Float64[]
    for row in eachrow(df)
        xv = row[1]
        yv = row[2]
        if ismissing(xv) || ismissing(yv) || xv == "" || yv == ""
            continue
        end
        push!(x, xv isa AbstractString ? parse(Float64, xv) : Float64(xv))
        push!(y, yv isa AbstractString ? parse(Float64, yv) : Float64(yv))
    end
    return x, y
end

function finite_xy(x::AbstractVector, y::AbstractVector)
    keep = isfinite.(Float64.(x)) .& isfinite.(Float64.(y))
    return Float64.(x[keep]), Float64.(y[keep])
end

function ensure_legacy_circuit_crop()
    img = FileIO.load(LEGACY_CIRCUIT_SOURCE)
    crop = img[61:420, 8:305]
    FileIO.save(LEGACY_CIRCUIT_CROP, crop)
    return LEGACY_CIRCUIT_CROP
end

function add_panel_label!(ax::Axis, label::String)
    GLMakie.text!(ax, 0.01, 0.99; text = label, space = :relative, align = (:left, :top), fontsize = 26, font = :bold, color = :black)
end

function add_l_scale_bar!(ax::Axis, x_min::Float64, x_max::Float64, y_min::Float64, y_max::Float64; time_len_s::Float64 = 20.0, volt_len_mv::Float64 = 50.0)
    x_span = x_max - x_min
    y_span = y_max - y_min
    x1 = x_max - 0.06 * x_span
    x0 = x1 - time_len_s
    y0 = y_min + 0.12 * y_span
    y1 = y0 + volt_len_mv
    GLMakie.lines!(ax, [x0, x0], [y0, y1]; color = :black, linewidth = 3)
    GLMakie.lines!(ax, [x0, x1], [y0, y0]; color = :black, linewidth = 3)
    GLMakie.text!(ax, x0 + 0.03 * x_span, (y0 + y1) / 2; text = "$(Int(round(volt_len_mv))) mV", rotation = pi / 2, align = (:center, :center), fontsize = 15, color = :black)
    GLMakie.text!(ax, (x0 + x1) / 2, y0 - 0.06 * y_span; text = "$(Int(round(time_len_s))) s", align = (:center, :center), fontsize = 15, color = :black)
end

function plot_circuit_panel!(ax::Axis)
    hidedecorations!(ax)
    hidespines!(ax)
    ax.aspect = DataAspect()
    img = rotr90(FileIO.load(ensure_legacy_circuit_crop()))
    GLMakie.image!(ax, img)
    ax.title = "Circuit Diagram"
end

function plot_four_trace_panel!(ax::Axis, t1, v1, t2, v2, t3, v3, t4, v4, label1, label2, label3, label4)
    offsets = [285.0, 190.0, 95.0, 0.0]
    time1 = Float64.(t1)
    time2 = Float64.(t2)
    time3 = Float64.(t3)
    time4 = Float64.(t4)
    vals1 = Float64.(v1)
    vals2 = Float64.(v2)
    vals3 = Float64.(v3)
    vals4 = Float64.(v4)
    GLMakie.lines!(ax, time1, vals1 .+ offsets[1]; color = RGBf(0.35, 0.35, 0.35), linewidth = 1.4)
    GLMakie.lines!(ax, time2, vals2 .+ offsets[2]; color = RGBf(0.1, 0.1, 0.1), linewidth = 1.4)
    GLMakie.lines!(ax, time3, vals3 .+ offsets[3]; color = MODE_COLORS["presynaptic"], linewidth = 1.4)
    GLMakie.lines!(ax, time4, vals4 .+ offsets[4]; color = MODE_COLORS["postsynaptic"], linewidth = 1.4)
    x_min = minimum((first(time1), first(time2), first(time3), first(time4)))
    x_max = maximum((last(time1), last(time2), last(time3), last(time4)))
    y_min = minimum(vals4)
    y_max = maximum(vals1 .+ offsets[1])
    label_x = x_min + 0.015 * (x_max - x_min)
    GLMakie.text!(ax, label_x, offsets[1] + 38; text = label1, align = (:left, :center), fontsize = 16, color = RGBf(0.35, 0.35, 0.35))
    GLMakie.text!(ax, label_x, offsets[2] + 38; text = label2, align = (:left, :center), fontsize = 16, color = RGBf(0.1, 0.1, 0.1))
    GLMakie.text!(ax, label_x, offsets[3] + 38; text = label3, align = (:left, :center), fontsize = 16, color = MODE_COLORS["presynaptic"])
    GLMakie.text!(ax, label_x, offsets[4] + 38; text = label4, align = (:left, :center), fontsize = 16, color = MODE_COLORS["postsynaptic"])
    ax.title = "Representative Traces"
    ax.ylabel = "V (mV)"
    ax.xticksvisible = false
    ax.yticksvisible = false
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xautolimitmargin = (0.0f0, 0.0f0)
    ax.yautolimitmargin = (0.0f0, 0.0f0)
    hidespines!(ax, :t, :r)
    GLMakie.xlims!(ax, x_min, x_max)
    GLMakie.ylims!(ax, y_min, y_max)
    add_l_scale_bar!(ax, x_min, x_max, y_min, y_max; time_len_s = 20.0, volt_len_mv = 50.0)
end

function plot_gain_panel!(ax::Axis, summary::DataFrame, raw::DataFrame, phase::String)
    ax.xlabel = "Gain"
    ax.ylabel = phase == "pre" ? "Freq on (Hz)" : "Freq off (Hz)"
    ax.title = phase == "pre" ? "During Drive" : "After Drive Removal"
    GLMakie.xlims!(ax, 0, 6)
    GLMakie.ylims!(ax, phase == "post" ? -0.005 : 0.0, 0.21)
    legend_entries = Any[]
    legend_labels = String[]
    for mode_name in ["presynaptic", "postsynaptic"]
        color = MODE_COLORS[mode_name]
        marker = MODE_MARKERS[mode_name]
        sub_summary = sort(summary[summary.mode .== mode_name, :], :control_gain)
        sub_raw = raw[(raw.mode .== mode_name) .& (raw.phase .== phase), :]
        bgx, bgy = finite_xy(sub_raw.control_gain, sub_raw.frequency_hz)
        ycol = phase == "pre" ? :mean_pre_frequency_hz : :mean_post_frequency_hz
        x, y = finite_xy(sub_summary.control_gain, sub_summary[!, ycol])
        !isempty(bgx) && GLMakie.scatter!(ax, bgx, bgy; color = (color, 0.24), marker = :circle, markersize = 7, strokewidth = 0)
        sc = GLMakie.scatter!(ax, x, y; color = color, marker = marker, markersize = 11, strokecolor = :black, strokewidth = 0.9)
        push!(legend_entries, sc)
        push!(legend_labels, clean_name(mode_name))
    end
    axislegend(ax, legend_entries, legend_labels; position = :rb, framevisible = true, backgroundcolor = (:white, 0.9), labelsize = 14, patchsize = (18, 18))
end

function plot_compare_panel!(ax::Axis, model_points::DataFrame)
    bio_x, bio_y = load_biological_fig2c()
    bio_handle = GLMakie.scatter!(ax, bio_x, bio_y; color = (:black, 0.9), marker = :diamond, markersize = 11, strokewidth = 0)
    handles = Any[bio_handle]
    labels = String["biology"]
    for mode_name in ["presynaptic", "postsynaptic"]
        sub = model_points[model_points.mode .== mode_name, :]
        handle = GLMakie.scatter!(ax, sub.si1_frequency_hz, sub.burst_frequency_hz; color = (MODE_COLORS[mode_name], 0.84), marker = MODE_MARKERS[mode_name], markersize = 11, strokecolor = :black, strokewidth = 0.8)
        push!(handles, handle)
        push!(labels, mode_name)
    end
    ax.title = "Si1 Rate vs Burst Rate"
    ax.xlabel = "Si1 or V0 spike freq (Hz)"
    ax.ylabel = "Burst freq (Hz)"
    GLMakie.xlims!(ax, 0, 6.5)
    GLMakie.ylims!(ax, 0, 0.21)
    hidespines!(ax, :t, :r)
    axislegend(ax, handles, labels; position = :rb, framevisible = true, backgroundcolor = (:white, 0.9), labelsize = 14, patchsize = (18, 18))
end

function main()
    for path in [PRESYNAPTIC_POINTS_CSV, PRESYNAPTIC_SUMMARY_CSV, POSTSYNAPTIC_POINTS_CSV, POSTSYNAPTIC_SUMMARY_CSV, BIO_DRIVEN_TRACES_CSV, BIO_DRIVEN_POINTS_CSV]
        isfile(path) || error("Missing $(path). Run the generating scripts first.")
    end

    raw = vcat(CSV.read(PRESYNAPTIC_POINTS_CSV, DataFrame), CSV.read(POSTSYNAPTIC_POINTS_CSV, DataFrame))
    summary = vcat(CSV.read(PRESYNAPTIC_SUMMARY_CSV, DataFrame), CSV.read(POSTSYNAPTIC_SUMMARY_CSV, DataFrame))
    model_points = CSV.read(BIO_DRIVEN_POINTS_CSV, DataFrame)
    driven_traces = CSV.read(BIO_DRIVEN_TRACES_CSV, DataFrame)

    bio_t1, bio_v1 = load_two_column_txt(BIO_SI1_TXT)
    bio_t2, bio_v2 = load_two_column_txt(BIO_SI2_TXT)
    bio_t1 .-= first(bio_t1)
    bio_t2 .-= first(bio_t2)

    pre_driven = driven_traces[driven_traces.mode .== "presynaptic", :]
    post_driven = driven_traces[driven_traces.mode .== "postsynaptic", :]

    fig = GLMakie.Figure(size = (1180, 1320), backgroundcolor = :white, fontsize = 20)
    axA = Axis(fig[1, 1], backgroundcolor = RGBf(0.98, 0.98, 0.99))
    plot_circuit_panel!(axA)
    axB = Axis(fig[1, 2], backgroundcolor = RGBf(0.98, 0.98, 0.99))
    plot_compare_panel!(axB, model_points)
    axC = Axis(fig[2, 1], backgroundcolor = RGBf(0.98, 0.98, 0.99))
    plot_gain_panel!(axC, summary, raw, "pre")
    axD = Axis(fig[2, 2], backgroundcolor = RGBf(0.98, 0.98, 0.99))
    plot_gain_panel!(axD, summary, raw, "post")
    axE = Axis(fig[3, 1:2], backgroundcolor = RGBf(0.98, 0.98, 0.99))
    plot_four_trace_panel!(
        axE,
        bio_t1, bio_v1,
        bio_t2, bio_v2,
        Float64.(pre_driven.time_s), Float64.(pre_driven.V1),
        Float64.(post_driven.time_s), Float64.(post_driven.V1),
        "Si1 biology", "Si2 biology", "Si2 presynaptic", "Si2 postsynaptic",
    )

    add_panel_label!(axA, "A")
    add_panel_label!(axB, "B")
    add_panel_label!(axC, "C")
    add_panel_label!(axD, "D")
    add_panel_label!(axE, "E")

    rowsize!(fig.layout, 1, Fixed(235))
    rowsize!(fig.layout, 2, Fixed(335))
    rowsize!(fig.layout, 3, Fixed(255))
    GLMakie.rowgap!(fig.layout, 18)
    GLMakie.colgap!(fig.layout, 16)
    GLMakie.resize_to_layout!(fig)

    mkpath(OUTPUT_DIR)
    CairoMakie.save(FIGURE_PNG, fig)
    CairoMakie.save(PAPER_FIG5_JPG, fig)
    display(fig)
    println("Saved $(FIGURE_PNG)")
    println("Saved $(PAPER_FIG5_JPG)")
end

main()
