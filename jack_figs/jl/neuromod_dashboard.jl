using CSV
using DataFrames
using GLMakie

include(joinpath(@__DIR__, "legacy", "compare_synaptic_plasticity_models.jl"))

const DASHBOARD_OUTPUT_DIR = joinpath(@__DIR__, "outputs")
const DASHBOARD_POINTS_CSV = joinpath(DASHBOARD_OUTPUT_DIR, "dashboard_scan_points.csv")
const DASHBOARD_SUMMARY_CSV = joinpath(DASHBOARD_OUTPUT_DIR, "dashboard_scan_summary.csv")

const DASHBOARD_PARAMS = updated_params(
    default_params();
    direct_post_base_g = 0.00064,
)

const DASHBOARD_SPECS = [
    ModeSweepSpec(presynaptic, "Presynaptic", vcat([0.0], collect(range(0.3050847457627119, 6.0; length = 119))), "Gain"),
    ModeSweepSpec(postsynaptic, "Postsynaptic", vcat([0.0], collect(range(0.06060606060606061, 6.0; length = 99))), "Gain"),
]

const DASHBOARD_COLORS = Dict(
    "presynaptic" => RGBf(0.12, 0.42, 0.86),
    "postsynaptic" => RGBf(0.80, 0.22, 0.18),
)
const DASHBOARD_MARKERS = Dict(
    "presynaptic" => :circle,
    "postsynaptic" => :rect,
)

function dashboard_config()
    base = default_config()
    SweepConfig(
        base.max_time_s,
        5.0,
        base.spike_threshold_mv,
        base.spike_refractory_s,
        base.burst_factor,
        base.onset_timeout_s,
        base.transient_bursts,
        base.measured_cycles,
    )
end

function empty_trajectory()
    (
        time_s = [0.0],
        V0 = [NaN],
        V3 = [NaN],
        V4 = [NaN],
        onset_failed = false,
    )
end

function selection_summary_row(summary::DataFrame, mode_name::String, control_gain::Float64)
    rows = summary[(summary.mode .== mode_name) .& (summary.control_gain .== control_gain), :]
    isempty(rows) && return nothing
    rows[1, :]
end

function selected_trajectory(mode_name::String, control_gain::Float64)
    sim = run_simulation(DASHBOARD_PARAMS, control_gain, parse_mode_label(mode_name), dashboard_config())
    (
        time_s = sim.time_s,
        V0 = sim.V0,
        V3 = sim.V3,
        V4 = sim.V4,
        onset_failed = sim.onset_failed,
    )
end

function main()
    isfile(DASHBOARD_POINTS_CSV) || error("Missing $(DASHBOARD_POINTS_CSV). Run dashboard_scans.jl first.")
    isfile(DASHBOARD_SUMMARY_CSV) || error("Missing $(DASHBOARD_SUMMARY_CSV). Run dashboard_scans.jl first.")

    raw = CSV.read(DASHBOARD_POINTS_CSV, DataFrame)
    summary = CSV.read(DASHBOARD_SUMMARY_CSV, DataFrame)

    start_row = sort(summary, [:mode, :control_gain])[1, :]
    selection = GLMakie.Observable((String(start_row.mode), Float64(start_row.control_gain)))
    row_obs = GLMakie.Observable(start_row)
    traj_obs = GLMakie.Observable(selected_trajectory(selection[][1], selection[][2]))

    GLMakie.on(selection) do (mode_name, gain)
        row_obs[] = selection_summary_row(summary, mode_name, gain)
        traj_obs[] = selected_trajectory(mode_name, gain)
    end

    fig = GLMakie.Figure(size = (1600, 1300), backgroundcolor = :white, fontsize = 18)
    title_text = GLMakie.lift(selection, row_obs, traj_obs) do sel, row, traj
        if isnothing(row)
            return "Selected: $(sel[1]) gain=$(round(sel[2], digits = 4))"
        end
        "Selected: $(sel[1]) gain=$(round(sel[2], digits = 4)) | freq on=$(round(row.mean_pre_frequency_hz, digits = 4)) Hz | freq off=$(round(row.mean_post_frequency_hz, digits = 4)) Hz" *
        (traj.onset_failed ? " | onset timeout" : "")
    end
    GLMakie.Label(fig[0, 1:2], title_text; tellwidth = false, fontsize = 24, justification = :left)

    pick_targets = IdDict{Any, Tuple{String, DataFrame, Symbol}}()
    high_x = Dict{Tuple{String, Symbol}, GLMakie.Observable{Vector{Float64}}}()
    high_y = Dict{Tuple{String, Symbol}, GLMakie.Observable{Vector{Float64}}}()

    panel_defs = [
        (:mean_pre_frequency_hz, "During Drive", "Freq on (Hz)"),
        (:mean_post_frequency_hz, "After Drive Removal", "Freq off (Hz)"),
    ]

    for (row_idx, spec) in enumerate(DASHBOARD_SPECS)
        mode_name = mode_label(spec.mode)
        sub_summary = sort(summary[summary.mode .== mode_name, :], :control_gain)
        sub_raw_pre = raw[(raw.mode .== mode_name) .& (raw.phase .== "pre"), :]
        sub_raw_post = raw[(raw.mode .== mode_name) .& (raw.phase .== "post"), :]

        for (col_idx, (ycol, title, ylabel)) in enumerate(panel_defs)
            ax = GLMakie.Axis(fig[row_idx, col_idx], title = "$(spec.label): $(title)", xlabel = "Gain", ylabel = ylabel)
            sub_raw = ycol == :mean_pre_frequency_hz ? sub_raw_pre : sub_raw_post
            !isempty(sub_raw) && GLMakie.scatter!(ax, Float64.(sub_raw.control_gain), Float64.(sub_raw.frequency_hz); color = (DASHBOARD_COLORS[mode_name], 0.24), markersize = 7)
            sc = GLMakie.scatter!(ax, Float64.(sub_summary.control_gain), Float64.(sub_summary[!, ycol]); color = DASHBOARD_COLORS[mode_name], marker = DASHBOARD_MARKERS[mode_name], markersize = 12, strokecolor = :black, strokewidth = 0.9)
            pick_targets[sc] = (mode_name, sub_summary, ycol)
            key = (mode_name, ycol)
            high_x[key] = GLMakie.Observable([Float64(sub_summary[1, :control_gain])])
            high_y[key] = GLMakie.Observable([Float64(sub_summary[1, ycol])])
            GLMakie.scatter!(ax, high_x[key], high_y[key]; color = :gold, markersize = 18, strokecolor = :black, strokewidth = 1.5, inspectable = false)
        end
    end

    axV0 = GLMakie.Axis(fig[3, 1:2], title = "V0", xlabel = "Time (s)", ylabel = "mV")
    axV4 = GLMakie.Axis(fig[4, 1:2], title = "V4", xlabel = "Time (s)", ylabel = "mV")
    axV3 = GLMakie.Axis(fig[5, 1:2], title = "V3", xlabel = "Time (s)", ylabel = "mV")

    line_t = GLMakie.lift(traj_obs) do tr
        tr.time_s
    end
    line_v0 = GLMakie.lift(traj_obs) do tr
        tr.V0
    end
    line_v4 = GLMakie.lift(traj_obs) do tr
        tr.V4
    end
    line_v3 = GLMakie.lift(traj_obs) do tr
        tr.V3
    end

    GLMakie.lines!(axV0, line_t, line_v0; color = :black, linewidth = 1.6)
    GLMakie.lines!(axV4, line_t, line_v4; color = RGBf(0.2, 0.2, 0.8), linewidth = 1.6)
    GLMakie.lines!(axV3, line_t, line_v3; color = RGBf(0.8, 0.2, 0.2), linewidth = 1.6)

    function refresh_axes!()
        tr = traj_obs[]
        tmin = minimum(tr.time_s)
        tmax = maximum(tr.time_s)
        tmax > tmin || (tmax = tmin + 1.0)
        for (ax, vals) in [(axV0, tr.V0), (axV4, tr.V4), (axV3, tr.V3)]
            vmin = minimum(vals)
            vmax = maximum(vals)
            if !isfinite(vmin) || !isfinite(vmax) || vmax <= vmin
                vmin, vmax = -80.0, 40.0
            else
                pad = 0.08 * (vmax - vmin)
                vmin -= pad
                vmax += pad
            end
            GLMakie.xlims!(ax, tmin, tmax)
            GLMakie.ylims!(ax, vmin, vmax)
        end
    end

    GLMakie.on(traj_obs) do _
        refresh_axes!()
    end
    refresh_axes!()

    GLMakie.on(GLMakie.events(fig).mousebutton) do event
        event.button == Mouse.left || return
        event.action == Mouse.press || return
        plt, idx = GLMakie.pick(fig)
        idx < 1 && return
        haskey(pick_targets, plt) || return
        mode_name, sub_summary, ycol = pick_targets[plt]
        gain = Float64(sub_summary[idx, :control_gain])
        selection[] = (mode_name, gain)
        for (m, col) in keys(high_x)
            if m == mode_name && col == ycol
                high_x[(m, col)][] = [gain]
                high_y[(m, col)][] = [Float64(sub_summary[idx, col])]
            else
                other = sort(summary[summary.mode .== m, :], :control_gain)
                if !isempty(other)
                    sel = selection_summary_row(summary, selection[][1], selection[][2])
                    if m == selection[][1]
                        high_x[(m, col)][] = [selection[][2]]
                        high_y[(m, col)][] = [Float64(sel[col])]
                    end
                end
            end
        end
    end

    GLMakie.display(fig)
end

main()
