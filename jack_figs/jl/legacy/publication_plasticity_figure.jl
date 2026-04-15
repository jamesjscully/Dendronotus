using GLMakie
using CSV
using DataFrames

include(joinpath(@__DIR__, "publication_plasticity_scan.jl"))

const PUB_PNG = joinpath(OUTPUT_DIR, "plasticity_publication_figure.png")
const AKIRA_FIG2C_TSV = joinpath(dirname(@__DIR__), "akira", "xlsx_cache", "Fig_2C.tsv")

function publication_modes()
    [
        ("presynaptic", "Presynaptic", RGBf(0.12, 0.42, 0.86), :circle),
        ("postsynaptic", "Postsynaptic", RGBf(0.80, 0.22, 0.18), :rect),
    ]
end

function publication_summary()
    summary = load_publication_summary()
    subset = summary[in.(summary[!, :mode], Ref(["presynaptic", "postsynaptic"])), :]
    isempty(subset) && error("No presynaptic/postsynaptic rows found in summary CSV")
    return subset
end

function publication_raw()
    raw = load_publication_raw()
    subset = raw[in.(raw[!, :mode], Ref(["presynaptic", "postsynaptic"])), :]
    isempty(subset) && error("No presynaptic/postsynaptic rows found in publication raw CSV")
    return subset
end

function biological_fig2c()
    isfile(AKIRA_FIG2C_TSV) || error("Missing Fig 2C TSV at $(AKIRA_FIG2C_TSV)")
    df = DataFrame(CSV.File(AKIRA_FIG2C_TSV; delim = '\t', normalizenames = false))
    x = Float64[]
    y = Float64[]
    for row in eachrow(df)
        x_raw = row[1]
        y_raw = row[2]
        if ismissing(x_raw) || ismissing(y_raw) || x_raw == "" || y_raw == ""
            continue
        end
        try
            push!(x, x_raw isa AbstractString ? parse(Float64, x_raw) : Float64(x_raw))
            push!(y, y_raw isa AbstractString ? parse(Float64, y_raw) : Float64(y_raw))
        catch
        end
    end
    return x, y
end

function finite_xy(x::AbstractVector, y::AbstractVector)
    keep = isfinite.(Float64.(x)) .& isfinite.(Float64.(y))
    return Float64.(x[keep]), Float64.(y[keep])
end

function compact_publication_figure(summary, raw)
    GLMakie.activate!()

    fig = GLMakie.Figure(size = (1700, 460), fontsize = 20, backgroundcolor = :white)
    ax_on = GLMakie.Axis(
        fig[1, 1];
        title = "A. During Drive",
        xlabel = "Gain",
        ylabel = "Freq on (Hz)",
        limits = (0, 6, 0, 0.21),
        xticks = 0:1:6,
        topspinevisible = true,
        rightspinevisible = true,
    )
    ax_off = GLMakie.Axis(
        fig[1, 2];
        title = "B. After Drive Removal",
        xlabel = "Gain",
        ylabel = "Freq off (Hz)",
        limits = (0, 6, 0, 0.21),
        xticks = 0:1:6,
        topspinevisible = true,
        rightspinevisible = true,
    )
    ax_compare = GLMakie.Axis(
        fig[1, 3];
        title = "C. Si1/V0 Rate vs Burst Rate",
        xlabel = "Si1 or V0 spike freq (Hz)",
        ylabel = "Burst freq (Hz)",
        limits = (0, 6.5, 0, 0.21),
        topspinevisible = true,
        rightspinevisible = true,
    )

    legend_entries = Any[]
    legend_labels = String[]

    bio_x, bio_y = biological_fig2c()
    bio_scatter = GLMakie.scatter!(
        ax_compare,
        bio_x,
        bio_y;
        color = (:black, 0.9),
        marker = :diamond,
        markersize = 11,
        strokewidth = 0,
    )
    push!(legend_entries, bio_scatter)
    push!(legend_labels, "Biology (Fig. 2C)")

    for (mode_name, label, color, marker) in publication_modes()
        subset = sort(summary[summary[!, :mode] .== mode_name, :], :control_gain)
        raw_subset = raw[raw[!, :mode] .== mode_name, :]
        raw_pre = raw_subset[raw_subset[!, :phase] .== "pre", :]
        raw_post = raw_subset[raw_subset[!, :phase] .== "post", :]

        x_on_bg, y_on_bg = finite_xy(raw_pre[!, :control_gain], raw_pre[!, :frequency_hz])
        x_off_bg, y_off_bg = finite_xy(raw_post[!, :control_gain], raw_post[!, :frequency_hz])
        x_on, y_on = finite_xy(subset[!, :control_gain], subset[!, :mean_pre_frequency_hz])
        x_off, y_off = finite_xy(subset[!, :control_gain], subset[!, :mean_post_frequency_hz])
        x_cmp_bg, y_cmp_bg = finite_xy(raw_pre[!, :driver_frequency_hz], raw_pre[!, :frequency_hz])
        x_cmp, y_cmp = finite_xy(subset[!, :mean_pre_driver_frequency_hz], subset[!, :mean_pre_frequency_hz])

        !isempty(x_on) || continue

        if !isempty(x_on_bg)
            GLMakie.scatter!(
                ax_on,
                x_on_bg,
                y_on_bg;
                color = (color, 0.12),
                marker = :circle,
                markersize = 7,
                strokewidth = 0,
            )
        end
        scatter_on = GLMakie.scatter!(
            ax_on,
            x_on,
            y_on;
            color = color,
            marker = marker,
            markersize = 17,
            strokecolor = :black,
            strokewidth = 1.0,
        )
        if !isempty(x_cmp_bg)
            GLMakie.scatter!(
                ax_compare,
                x_cmp_bg,
                y_cmp_bg;
                color = (color, 0.12),
                marker = :circle,
                markersize = 7,
                strokewidth = 0,
            )
        end
        GLMakie.scatter!(
            ax_compare,
            x_cmp,
            y_cmp;
            color = color,
            marker = marker,
            markersize = 17,
            strokecolor = :black,
            strokewidth = 1.0,
        )

        if !isempty(x_off_bg)
            GLMakie.scatter!(
                ax_off,
                x_off_bg,
                y_off_bg;
                color = (color, 0.12),
                marker = :circle,
                markersize = 7,
                strokewidth = 0,
            )
        end
        GLMakie.scatter!(
            ax_off,
            x_off,
            y_off;
            color = color,
            marker = marker,
            markersize = 17,
            strokecolor = :black,
            strokewidth = 1.0,
        )

        push!(legend_entries, scatter_on)
        push!(legend_labels, label)
    end

    GLMakie.Legend(fig[0, 1:3], legend_entries, legend_labels; orientation = :horizontal, tellwidth = false, framevisible = false)

    GLMakie.rowgap!(fig.layout, 8)
    GLMakie.colgap!(fig.layout, 18)
    return fig
end

function main()
    mkpath(OUTPUT_DIR)
    if !isfile(PUB_SCAN_SUMMARY_CSV) || !isfile(PUB_SCAN_RAW_CSV)
        raw, summary = generate_publication_data()
        save_publication_data(raw, summary)
    end
    summary = publication_summary()
    raw = publication_raw()
    fig = compact_publication_figure(summary, raw)
    screen = GLMakie.Screen()
    GLMakie.display(screen, fig)
    GLMakie.save(PUB_PNG, fig)
    println("Saved:")
    println("  $(PUB_PNG)")
    wait(screen.scene)
end

const PUBFIG_NOAUTORUN = get(ENV, "PLASTICITY_PUBFIG_NOAUTORUN", "0") == "1"
const PUBFIG_SCRIPT_FILE = abspath(@__FILE__)

if !PUBFIG_NOAUTORUN && (abspath(PROGRAM_FILE) == PUBFIG_SCRIPT_FILE || isinteractive())
    main()
elseif isinteractive()
    @info "Publication figure loaded with autorun disabled. Run `main()` to render it."
end
