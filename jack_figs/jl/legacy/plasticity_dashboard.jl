using GLMakie
using Statistics

include(joinpath(@__DIR__, "compare_synaptic_plasticity_models.jl"))

const DASHBOARD_FIGURE = joinpath(OUTPUT_DIR, "plasticity_dashboard_snapshot.png")

function mode_colors()
    Dict(
        "direct_gmax" => RGBf(0.0, 0.45, 0.05),
        "presynaptic" => RGBf(0.12, 0.49, 0.89),
        "postsynaptic" => RGBf(0.83, 0.18, 0.18),
    )
end

function summary_or_generate(config::SweepConfig, params::ModelParams, specs::Vector{ModeSweepSpec})
    if isfile(SUMMARY_CSV)
        return load_summary()
    end

    raw_points, summary = generate_data(config, params, specs)
    save_data(raw_points, summary)
    return summary
end

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

function measured_window(bursts::BurstDetection, config::SweepConfig)
    required_starts = config.transient_bursts + config.measured_cycles + 1
    complete_starts = complete_burst_start_times(bursts)
    length(complete_starts) < required_starts && return (NaN, NaN)

    start_idx = config.transient_bursts + 1
    stop_idx = config.transient_bursts + config.measured_cycles + 1
    start_t = complete_starts[start_idx]
    stop_t = complete_starts[stop_idx]
    return start_t, stop_t
end

function summarize_selected_point(summary::DataFrame, mode_name::String, control_gain::Float64)
    rows = summary[(summary.mode .== mode_name) .& (summary.control_gain .== control_gain), :]
    isempty(rows) && return nothing
    return rows[1, :]
end

function selected_trajectory(mode_name::String, control_gain::Float64, params::ModelParams, config::SweepConfig)
    mode = parse_mode_label(mode_name)
    sim = run_simulation(params, control_gain, mode, config)
    traces = effective_excitation_trace(sim.sm, sim.s4, control_gain, mode, sim.params)
    bursts = detect_bursts_from_spikes(sim.spike_times_s, config)
    measured_start, measured_stop = measured_window(bursts, config)
    return (
        time_s = sim.time_s,
        V0 = sim.V0,
        V1 = sim.V1,
        V2 = sim.V2,
        V3 = sim.V3,
        V4 = sim.V4,
        sm = sim.sm,
        s4 = sim.s4,
        excitation = traces.strength,
        gain = traces.gain,
        spike_times_s = sim.spike_times_s,
        burst_starts_s = bursts.start_times_s,
        complete_burst_starts_s = complete_burst_start_times(bursts),
        incomplete_burst_starts_s = incomplete_burst_start_times(bursts),
        burst_ends_s = bursts.end_times_s,
        burst_factor = bursts.burst_factor,
        burst_counts = bursts.spike_counts,
        complete_burst_count = count(identity, bursts.is_complete),
        incomplete_burst_count = count(!, bursts.is_complete),
        measured_start_s = sim.measured_start_s,
        measured_stop_s = sim.measured_stop_s,
        stimulus_off_s = sim.stimulus_off_s,
        withdrawal_end_s = sim.withdrawal_end_s,
        onset_failed = sim.onset_failed,
        loaded = true,
    )
end

function find_spec(specs::Vector{ModeSweepSpec}, mode_name::String)
    idx = findfirst(spec -> mode_label(spec.mode) == mode_name, specs)
    idx === nothing && error("No spec found for mode $(mode_name)")
    return specs[idx]
end

function trajectory_cache_key(mode_name::String, control_gain::Float64)
    return (mode_name, round(control_gain; digits = 10))
end

function empty_trajectory()
    return (
        time_s = [0.0],
        V0 = [NaN],
        V1 = [NaN],
        V2 = [NaN],
        V3 = [NaN],
        V4 = [NaN],
        sm = [NaN],
        s4 = [NaN],
        excitation = [NaN],
        gain = [NaN],
        spike_times_s = Float64[],
        burst_starts_s = Float64[],
        complete_burst_starts_s = Float64[],
        incomplete_burst_starts_s = Float64[],
        burst_ends_s = Float64[],
        burst_factor = NaN,
        burst_counts = Int[],
        complete_burst_count = 0,
        incomplete_burst_count = 0,
        measured_start_s = NaN,
        measured_stop_s = NaN,
        stimulus_off_s = NaN,
        withdrawal_end_s = NaN,
        onset_failed = false,
        loaded = false,
    )
end

function offscreen_vlines(traj, values)
    isempty(values) ? [first(traj.time_s) - 1.0] : values
end

function build_dashboard(; save_snapshot::Bool = false, initial_mode_name::Union{Nothing, String} = nothing, initial_control_gain::Union{Nothing, Float64} = nothing)
    GLMakie.activate!()

    config = dashboard_config()
    params = default_params()
    specs = default_sweep_specs()
    summary = summary_or_generate(config, params, specs)
    isempty(summary) && error("No valid burst summaries found with the current detector. Refine the burst heuristic before launching the dashboard.")
    colors = mode_colors()
    spec_map = Dict(mode_label(spec.mode) => spec for spec in specs)

    fallback_row = sort(summary, [:mode, :control_gain])[1, :]
    start_mode = something(initial_mode_name, String(fallback_row.mode))
    start_gain = if isnothing(initial_control_gain)
        if start_mode == String(fallback_row.mode)
            Float64(fallback_row.control_gain)
        else
            spec = spec_map[start_mode]
            first(spec.control_values)
        end
    else
        initial_control_gain
    end
    selection = Observable((mode = start_mode, control_gain = Float64(start_gain)))
    selected_row = Observable(summarize_selected_point(summary, selection[].mode, selection[].control_gain))
    trajectory_cache = Dict{Tuple{String, Float64}, NamedTuple}()

    function fetch_trajectory(mode_name::String, control_gain::Float64)
        key = trajectory_cache_key(mode_name, control_gain)
        get!(trajectory_cache, key) do
            selected_trajectory(mode_name, control_gain, params, config)
        end
    end

    trajectory = Observable(empty_trajectory())

    on(selection) do sel
        selected_row[] = summarize_selected_point(summary, sel.mode, sel.control_gain)
        trajectory[] = fetch_trajectory(sel.mode, sel.control_gain)
    end

    if !isnothing(initial_mode_name) || !isnothing(initial_control_gain)
        trajectory[] = fetch_trajectory(selection[].mode, selection[].control_gain)
    end

    fig = GLMakie.Figure(size = (1850, 1900), fontsize = 18)
    title_text = lift(selection, selected_row, trajectory) do sel, row, traj
        spec = find_spec(specs, sel.mode)
        row_text = if isnothing(row)
            " | no summary point under current strict criterion"
        else
            string(
                " | pre frequency = ", round(row.mean_pre_burst_frequency_hz; digits = 4), " Hz",
                " | post bursts = ", round(row.post_withdrawal_burst_count; digits = 3),
                " | post/pre = ", round(row.post_pre_frequency_ratio; digits = 4),
            )
        end
        string(
            "Selected: ", spec.label,
            " | control = ", round(sel.control_gain; digits = 4),
            row_text,
            traj.loaded && traj.onset_failed ? " | onset timeout" : "",
            " | complete/incomplete bursts = ", traj.loaded ? string(traj.complete_burst_count, "/", traj.incomplete_burst_count) : "not loaded",
            " | burst factor = ", isnan(traj.burst_factor) ? "n/a" : string(round(traj.burst_factor; digits = 3)),
        )
    end
    GLMakie.Label(fig[0, 1:3], title_text; tellwidth = false, fontsize = 24, justification = :left)

    panel_defs = [
        (xcol = :control_gain, ycol = :mean_pre_burst_frequency_hz, xlabel = :control, ylabel = "Mean Pre-Withdrawal Burst Frequency (Hz)", suffix = "Pre Frequency vs Control"),
        (xcol = :control_gain, ycol = :post_withdrawal_burst_count, xlabel = :control, ylabel = "Post-Withdrawal Complete Burst Count", suffix = "Post Bursts vs Control"),
        (xcol = :control_gain, ycol = :post_pre_frequency_ratio, xlabel = :control, ylabel = "Post/Pre Burst Frequency Ratio", suffix = "Post/Pre Ratio vs Control"),
    ]

    pick_targets = IdDict{Any, Tuple{String, DataFrame, Symbol, Symbol}}()
    highlight_x = Dict{Tuple{String, Symbol, Symbol}, Observable{Vector{Float64}}}()
    highlight_y = Dict{Tuple{String, Symbol, Symbol}, Observable{Vector{Float64}}}()

    for (row_idx, spec) in enumerate(specs)
        mode_name = mode_label(spec.mode)
        subset = sort(summary[summary.mode .== mode_name, :], :control_gain)

        for (col_idx, panel) in enumerate(panel_defs)
            ax = GLMakie.Axis(
                fig[row_idx, col_idx];
                title = "$(spec.label): $(panel.suffix)",
                xlabel = spec.control_axis_label,
                ylabel = panel.ylabel,
            )

            xs = Float64.(subset[!, panel.xcol])
            ys = Float64.(subset[!, panel.ycol])
            isempty(xs) && continue
            scatter_plot = GLMakie.scatter!(
                ax,
                xs,
                ys;
                color = colors[mode_name],
                markersize = 14,
                strokecolor = :black,
                strokewidth = 1.0,
            )
            pick_targets[scatter_plot] = (mode_name, subset, panel.xcol, panel.ycol)

            key = (mode_name, panel.xcol, panel.ycol)
            highlight_x[key] = Observable([xs[1]])
            highlight_y[key] = Observable([ys[1]])
            GLMakie.scatter!(
                ax,
                highlight_x[key],
                highlight_y[key];
                color = :gold,
                markersize = 22,
                strokecolor = :black,
                strokewidth = 2.0,
                inspectable = false,
            )
        end
    end

    voltage_ax0 = GLMakie.Axis(
        fig[4, 1:3];
        title = "Selected Trajectory: V0",
        xlabel = "Time (s)",
        ylabel = "Voltage (mV)",
    )
    voltage_ax4 = GLMakie.Axis(
        fig[5, 1:3];
        title = "Selected Trajectory: V4",
        xlabel = "Time (s)",
        ylabel = "Voltage (mV)",
    )
    voltage_ax3 = GLMakie.Axis(
        fig[6, 1:3];
        title = "Selected Trajectory: V3",
        xlabel = "Time (s)",
        ylabel = "Voltage (mV)",
    )

    time_obs = lift(traj -> traj.time_s, trajectory)
    V0_obs = lift(traj -> traj.V0, trajectory)
    V3_obs = lift(traj -> traj.V3, trajectory)
    V4_obs = lift(traj -> traj.V4, trajectory)
    spike_times_obs = lift(traj -> offscreen_vlines(traj, traj.spike_times_s), trajectory)
    complete_burst_starts_obs = lift(traj -> offscreen_vlines(traj, traj.complete_burst_starts_s), trajectory)
    incomplete_burst_starts_obs = lift(traj -> offscreen_vlines(traj, traj.incomplete_burst_starts_s), trajectory)
    measured_start_obs = lift(traj -> offscreen_vlines(traj, isnan(traj.measured_start_s) ? Float64[] : [traj.measured_start_s]), trajectory)
    measured_stop_obs = lift(traj -> offscreen_vlines(traj, isnan(traj.measured_stop_s) ? Float64[] : [traj.measured_stop_s]), trajectory)
    stimulus_off_obs = lift(traj -> offscreen_vlines(traj, isnan(traj.stimulus_off_s) ? Float64[] : [traj.stimulus_off_s]), trajectory)
    withdrawal_end_obs = lift(traj -> offscreen_vlines(traj, isnan(traj.withdrawal_end_s) ? Float64[] : [traj.withdrawal_end_s]), trajectory)

    for (ax, signal_obs, color, label) in (
        (voltage_ax0, V0_obs, :black, "V0"),
        (voltage_ax4, V4_obs, :firebrick, "V4"),
        (voltage_ax3, V3_obs, :royalblue, "V3"),
    )
        GLMakie.lines!(ax, time_obs, signal_obs; color = color, linewidth = 1.5, label = label)
        GLMakie.vlines!(ax, spike_times_obs; color = (:gray, 0.18), linewidth = 0.8)
        GLMakie.vlines!(ax, complete_burst_starts_obs; color = (:black, 0.65), linewidth = 1.6)
        GLMakie.vlines!(ax, incomplete_burst_starts_obs; color = (:orange, 0.9), linewidth = 2.0)
        GLMakie.vlines!(ax, measured_start_obs; color = (:black, 0.8), linewidth = 2.0, linestyle = :dash)
        GLMakie.vlines!(ax, measured_stop_obs; color = (:black, 0.8), linewidth = 2.0, linestyle = :dash)
        GLMakie.vlines!(ax, stimulus_off_obs; color = (:purple, 0.9), linewidth = 2.2, linestyle = :dashdot)
        GLMakie.vlines!(ax, withdrawal_end_obs; color = (:purple, 0.7), linewidth = 1.8, linestyle = :dot)
        GLMakie.axislegend(ax; position = :rb)
    end

    function focus_measured_window!(traj)
        if !traj.loaded
            return
        elseif isnan(traj.measured_start_s) || isnan(traj.measured_stop_s)
            t_left = first(traj.time_s)
            t_right = last(traj.time_s)
        else
            cycle_span = max((traj.measured_stop_s - traj.measured_start_s) / max(config.measured_cycles, 1), 0.5)
            t_left = max(first(traj.time_s), traj.measured_start_s - 1.5 * cycle_span)
            t_right = min(last(traj.time_s), traj.measured_stop_s + 0.75 * cycle_span)
        end

        if t_right <= t_left
            t_right = t_left + 1.0
        end

        for (ax, signal) in (
            (voltage_ax0, traj.V0),
            (voltage_ax4, traj.V4),
            (voltage_ax3, traj.V3),
        )
            GLMakie.xlims!(ax, t_left, t_right)

            voltage_idx = findall(t -> t >= t_left && t <= t_right, traj.time_s)
            if isempty(voltage_idx)
                voltage_idx = 1:length(traj.time_s)
            end

            v_window = filter(isfinite, signal[voltage_idx])
            if !isempty(v_window)
                v_min = minimum(v_window)
                v_max = maximum(v_window)
                v_pad = max(2.0, 0.08 * (v_max - v_min + eps()))
                GLMakie.ylims!(ax, v_min - v_pad, v_max + v_pad)
            end
        end
    end
    focus_measured_window!(trajectory[])
    on(trajectory) do traj
        focus_measured_window!(traj)
    end

    function sync_highlights!(sel)
        row = summarize_selected_point(summary, sel.mode, sel.control_gain)
        isnothing(row) && return
        for panel in panel_defs
            key = (sel.mode, panel.xcol, panel.ycol)
            haskey(highlight_x, key) || continue
            highlight_x[key][] = [Float64(row[panel.xcol])]
            highlight_y[key][] = [Float64(row[panel.ycol])]
        end
    end
    sync_highlights!(selection[])
    on(selection) do sel
        sync_highlights!(sel)
    end

    on(GLMakie.events(fig).mousebutton) do event
        if event.button != GLMakie.Mouse.left || event.action != GLMakie.Mouse.press
            return
        end

        picked = GLMakie.pick(fig)
        picked === nothing && return
        plt, idx = picked
        haskey(pick_targets, plt) || return
        idx isa Integer || return
        idx < 1 && return

        mode_name, subset, _, _ = pick_targets[plt]
        idx > nrow(subset) && return
        row = subset[idx, :]
        selection[] = (mode = mode_name, control_gain = Float64(row.control_gain))
    end

    if save_snapshot
        GLMakie.save(DASHBOARD_FIGURE, fig)
    else
        screen = GLMakie.Screen()
        GLMakie.display(screen, fig)
        wait(screen.scene)
    end

    return fig
end

function main(args = ARGS)
    save_snapshot = any(arg -> lowercase(arg) == "snapshot", args)
    filtered_args = [arg for arg in args if lowercase(arg) != "snapshot"]
    initial_mode_name = nothing
    initial_control_gain = nothing

    if length(filtered_args) >= 1
        initial_mode_name = filtered_args[1]
    end
    if length(filtered_args) >= 2
        initial_control_gain = parse(Float64, filtered_args[2])
    end

    build_dashboard(; save_snapshot = save_snapshot, initial_mode_name = initial_mode_name, initial_control_gain = initial_control_gain)
end

launch_dashboard(; save_snapshot::Bool = false) = build_dashboard(; save_snapshot = save_snapshot)

const DASHBOARD_NOAUTORUN = get(ENV, "PLASTICITY_DASHBOARD_NOAUTORUN", "0") == "1"
const DASHBOARD_SCRIPT_FILE = abspath(@__FILE__)

if !DASHBOARD_NOAUTORUN && (abspath(PROGRAM_FILE) == DASHBOARD_SCRIPT_FILE || isinteractive())
    main()
elseif isinteractive()
    @info "Plasticity dashboard loaded with autorun disabled. Run `launch_dashboard()` to open it."
end
