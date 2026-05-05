using CairoMakie
using FileIO
using GLMakie
using Statistics

import GLMakie:
    Axis,
    DataAspect,
    Button,
    Checkbox,
    Figure,
    Fixed,
    GridLayout,
    Label,
    Menu,
    Mouse,
    Observable,
    Relative,
    RGBf,
    SliderGrid,
    colgap!,
    colsize!,
    display,
    events,
    hidedecorations!,
    hidespines!,
    image!,
    lift,
    lines!,
    pick,
    rowgap!,
    rowsize!,
    scatter!,
    set_close_to!,
    text!,
    xlims!,
    ylims!

if !isdefined(@__MODULE__, :dashboard_params)
    include(joinpath(@__DIR__, "neuromod_dashboard_numerics.jl"))
end

const DASH_COLORS = Dict(
    "presynaptic" => RGBf(0.12, 0.42, 0.86),
    "postsynaptic" => RGBf(0.80, 0.22, 0.18),
)
const DASH_MARKERS = Dict(
    "presynaptic" => :circle,
    "postsynaptic" => :rect,
)
const LITERATURE_GREY = RGBf(0.55, 0.55, 0.55)
const DASHBOARD_HANDLE = Ref{Any}(nothing)
const DASHBOARD_PARAM_KEYS = [
    :x_shift_si2,
    :x_shift_si3,
    :presyn_g0,
    :postsyn_g0,
    :alpha1,
    :beta1,
    :s0_floor,
    :alpham,
    :betam,
    :sm_floor,
    :ca_shift_si2,
    :ca_shift_si3,
    :presynaptic_base_g,
    :direct_post_base_g,
    :presyn_alphax,
    :postsyn_alphax,
    :presyn_betax,
    :postsyn_betax,
    :si3_exc_floor,
    :slow_inhib_g,
    :slow_inhib_alpha,
    :slow_inhib_beta,
    :si2_mutual_inhib_g,
    :si2_mutual_inhib_alpha,
    :si2_mutual_inhib_beta,
    :si_inhib_floor,
    :si3_mutual_inhib_g,
    :si3_mutual_inhib_alpha,
    :si3_mutual_inhib_beta,
    :t1_ms,
]

function built_in_param_defaults()
    return [
        CALIBRATED_X_SHIFT,
        CALIBRATED_X_SHIFT,
        PRESYNAPTIC_SI1_EXCITATORY_G,
        POSTSYNAPTIC_SI1_EXCITATORY_G,
        default_params().alpha1,
        default_params().beta1,
        default_params().s0_floor,
        default_params().alpham,
        default_params().betam,
        default_params().sm_floor,
        0.0,
        0.0,
        default_params().presynaptic_base_g,
        COMPENSATED_DIRECT_POST_BASE_G,
        PRESYNAPTIC_BASE_ALPHA,
        CITED_EXCITATORY_ALPHA,
        default_params().betax,
        default_params().betax,
        default_params().si3_exc_floor,
        default_params().g14,
        default_params().alphai,
        default_params().betai,
        default_params().g12,
        default_params().alpha2,
        default_params().beta2,
        default_params().si_inhib_floor,
        default_params().g34,
        default_params().alpha3,
        default_params().beta3,
        75000.0,
    ]
end

function load_param_defaults()
    defaults = Float64.(built_in_param_defaults())
    isfile(DASHBOARD_PARAM_DEFAULTS_CSV) || return defaults
    saved = CSV.read(DASHBOARD_PARAM_DEFAULTS_CSV, DataFrame)
    (:parameter in propertynames(saved) && :value in propertynames(saved)) || return defaults
    saved_values = Dict(Symbol(row.parameter) => Float64(row.value) for row in eachrow(saved))
    return [get(saved_values, key, defaults[i]) for (i, key) in enumerate(DASHBOARD_PARAM_KEYS)]
end

function save_param_defaults(values::AbstractVector{<:Real})
    mkpath(DASHBOARD_OUTPUT_DIR)
    df = DataFrame(parameter = String.(DASHBOARD_PARAM_KEYS), value = Float64.(values))
    CSV.write(DASHBOARD_PARAM_DEFAULTS_CSV, df)
    return DASHBOARD_PARAM_DEFAULTS_CSV
end

function trace_offsets()
    return Dict(
        "si1" => 360.0,
        "si2" => 265.0,
        "presynaptic" => 145.0,
        "presynaptic_si3" => 70.0,
        "postsynaptic" => -25.0,
        "postsynaptic_si3" => -100.0,
    )
end

shared_neuromod_offset(offsets) = offsets["si2"] - 58.0

function neuromod_trace_vectors(time_s, sm, base_offset::Float64)
    t, y = matched_trace_vectors(time_s, sm, 0.0)
    keep = isfinite.(y)
    if isempty(y) || !any(keep)
        return Float64[], Float64[]
    end
    return t[keep], base_offset .+ 34.0 .* clamp.(y[keep], 0.0, 1.5)
end

function add_panel_label!(ax::Axis, label::String)
    text!(ax, 0.01, 0.99; text = label, space = :relative, align = (:left, :top), fontsize = 24, font = :bold, color = :black)
end

function ensure_dashboard_circuit_crop()
    isfile(DASHBOARD_CIRCUIT_CROP) && return DASHBOARD_CIRCUIT_CROP
    img = FileIO.load(DASHBOARD_CIRCUIT_SOURCE)
    crop = img[61:420, 8:305]
    mkpath(DASHBOARD_OUTPUT_DIR)
    FileIO.save(DASHBOARD_CIRCUIT_CROP, crop)
    return DASHBOARD_CIRCUIT_CROP
end

function plot_circuit_panel!(ax::Axis)
    hidedecorations!(ax)
    hidespines!(ax)
    ax.aspect = DataAspect()
    image!(ax, rotr90(FileIO.load(ensure_dashboard_circuit_crop())))
    ax.title = "Circuit Diagram"
end

function model_point_vectors(points::DataFrame, mode_name::String)
    sub = points[points.mode .== mode_name, :]
    return finite_cols(sub, :si1_frequency_hz, :burst_frequency_hz)
end

function trend_line_vectors(x::Vector{Float64}, y::Vector{Float64})
    keep = isfinite.(x) .& isfinite.(y)
    x = x[keep]
    y = y[keep]
    length(x) < 2 && return Float64[], Float64[]
    x_min = minimum(x)
    x_max = maximum(x)
    x_max == x_min && return Float64[], Float64[]
    x_mean = mean(x)
    y_mean = mean(y)
    denom = sum((x .- x_mean) .^ 2)
    denom <= eps() && return Float64[], Float64[]
    slope = sum((x .- x_mean) .* (y .- y_mean)) / denom
    intercept = y_mean - slope * x_mean
    line_x = [x_min, x_max]
    return line_x, intercept .+ slope .* line_x
end

function matched_xy_vectors(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = min(length(x), length(y))
    n == 0 && return Float64[], Float64[]
    return Float64.(x[1:n]), Float64.(y[1:n])
end

function compare_panel_ylimit(model_points::DataFrame, bio_points::DataFrame)
    vals = Float64[]
    for mode_name in MODE_LABELS
        x, y = model_point_vectors(model_points, mode_name)
        append!(vals, y)
        _, ty = trend_line_vectors(x, y)
        append!(vals, ty)
    end
    bio_x, bio_y = finite_cols(bio_points, :si1_frequency_hz, :burst_frequency_hz)
    append!(vals, bio_y)
    _, bio_ty = trend_line_vectors(bio_x, bio_y)
    append!(vals, bio_ty)
    vals = vals[isfinite.(vals)]
    isempty(vals) && return 0.21
    return max(0.21, 1.12 * maximum(vals))
end

detector_kwargs(detector) = (
    spike_threshold_mv = Float64(detector.spike_threshold_mv),
    spike_refractory_s = Float64(detector.spike_refractory_s),
    burst_factor = Float64(detector.burst_factor),
)

function refresh_model_point_observables!(model_points_obs, point_obs)
    points = model_points_obs[]
    for mode_name in MODE_LABELS
        x, y = model_point_vectors(points, mode_name)
        point_obs[(mode_name, :x)][] = x
        point_obs[(mode_name, :y)][] = y
    end
end

function plot_compare_panel!(ax::Axis, model_points::DataFrame; detector = (spike_threshold_mv = -20.0, spike_refractory_s = 0.02, burst_factor = 3.0), time_scale::Float64 = REPRESENTATIVE_BIOLOGY_TIME_SCALE)
    lit_x, lit_y = load_biological_fig2c(; time_scale = time_scale)
    bio_points = load_biological_trace_points(; time_scale = time_scale, detector_kwargs(detector)...)
    bio_x, bio_y = finite_cols(bio_points, :si1_frequency_hz, :burst_frequency_hz)
    scatter!(ax, lit_x, lit_y; color = LITERATURE_GREY, marker = :circle, markersize = 10, strokecolor = :black, strokewidth = 0.8)
    scatter!(ax, bio_x, bio_y; color = :black, marker = :utriangle, markersize = 10, strokecolor = :white, strokewidth = 0.8)
    for mode_name in MODE_LABELS
        x, y = model_point_vectors(model_points, mode_name)
        scatter!(ax, x, y; color = (DASH_COLORS[mode_name], 0.84), marker = DASH_MARKERS[mode_name], markersize = 10, strokecolor = :black, strokewidth = 0.8)
    end
    lit_tx, lit_ty = trend_line_vectors(lit_x, lit_y)
    lines!(ax, lit_tx, lit_ty; color = LITERATURE_GREY, linewidth = 3.0)
    tx, ty = trend_line_vectors(bio_x, bio_y)
    lines!(ax, tx, ty; color = :black, linewidth = 4.0)
    for mode_name in MODE_LABELS
        x, y = model_point_vectors(model_points, mode_name)
        tx, ty = trend_line_vectors(x, y)
        lines!(ax, tx, ty; color = DASH_COLORS[mode_name], linewidth = 4.0)
    end
    ax.title = "Si1 Rate vs Burst Rate"
    ax.xlabel = "Si1 spike freq (Hz)"
    ax.ylabel = "Burst freq (Hz)"
    xlims!(ax, 0, 6.5)
    ylims!(ax, 0, compare_panel_ylimit(model_points, bio_points))
    hidespines!(ax, :t, :r)
end

function add_l_scale_bar!(ax::Axis, x_min::Float64, x_max::Float64, y_min::Float64, y_max::Float64; time_len_s::Float64 = 20.0, volt_len_mv::Float64 = 50.0)
    x_span = x_max - x_min
    y_span = y_max - y_min
    x1 = x_max - 0.06 * x_span
    x0 = x1 - time_len_s
    y0 = y_min + 0.12 * y_span
    y1 = y0 + volt_len_mv
    lines!(ax, [x0, x0], [y0, y1]; color = :black, linewidth = 3)
    lines!(ax, [x0, x1], [y0, y0]; color = :black, linewidth = 3)
    text!(ax, x0 + 0.03 * x_span, (y0 + y1) / 2; text = "$(Int(round(volt_len_mv))) mV", rotation = pi / 2, align = (:center, :center), fontsize = 14, color = :black)
    text!(ax, (x0 + x1) / 2, y0 - 0.06 * y_span; text = "$(Int(round(time_len_s))) s", align = (:center, :center), fontsize = 14, color = :black)
end

function burst_marker_vectors(time_s, voltage_mv, offset::Float64; spike_threshold_mv::Float64 = -20.0, spike_refractory_s::Float64 = 0.02, burst_factor::Float64 = 2.0)
    isempty(time_s) && return Float64[], Float64[]
    starts, _ = detect_bursts_from_trace_spikes(time_s, voltage_mv; spike_threshold_mv = spike_threshold_mv, spike_refractory_s = spike_refractory_s, burst_factor = burst_factor)
    isempty(starts) && return Float64[], Float64[]
    x = Float64.(time_s[starts])
    y = Float64.(voltage_mv[starts]) .+ offset
    return x, y
end

function burst_tick_vectors(time_s, voltage_mv, offset::Float64; spike_threshold_mv::Float64 = -20.0, spike_refractory_s::Float64 = 0.02, burst_factor::Float64 = 2.0)
    x, _ = burst_marker_vectors(time_s, voltage_mv, offset; spike_threshold_mv = spike_threshold_mv, spike_refractory_s = spike_refractory_s, burst_factor = burst_factor)
    isempty(x) && return Float64[], Float64[]
    xs = Float64[]
    ys = Float64[]
    for t in x
        append!(xs, [t, t, NaN])
        append!(ys, [offset - 48.0, offset + 52.0, NaN])
    end
    return xs, ys
end

function plot_burst_markers!(ax::Axis, time_s, voltage_mv, offset::Float64; color = :black, marker = :circle, strokecolor = :black, visible = true, detector = (spike_threshold_mv = -20.0, spike_refractory_s = 0.02, burst_factor = 3.0))
    tx, ty = burst_tick_vectors(time_s, voltage_mv, offset; detector_kwargs(detector)...)
    mx, my = burst_marker_vectors(time_s, voltage_mv, offset; detector_kwargs(detector)...)
    lines!(ax, tx, ty; color = (color, 0.28), linewidth = 1.0, visible = visible)
    scatter!(ax, mx, my; color = color, marker = marker, markersize = 10, strokecolor = strokecolor, strokewidth = 0.8, visible = visible)
end

function trace_panel_limits(bio_t1, bio_v1, bio_t2, bio_v2, traces::Dict{String, DataFrame}; source::Symbol = :representative, show_si3::Bool = false)
    xs = Float64[]
    ymins = Float64[]
    ymaxs = Float64[]
    offsets = trace_offsets()
    if source == :scan_protocol && !isempty(traces["presynaptic"]) && !isempty(traces["postsynaptic"])
        pre = traces["presynaptic"]
        append!(xs, [first(pre.time_s), last(pre.time_s)])
        pre_v0 = Float64.(pre.V0) .+ offsets["si1"]
        push!(ymins, minimum(pre_v0))
        push!(ymaxs, maximum(pre_v0))
    else
        append!(xs, [first(bio_t1), last(bio_t1), first(bio_t2), last(bio_t2)])
        push!(ymins, minimum(Float64.(bio_v1) .+ offsets["si1"]))
        push!(ymaxs, maximum(Float64.(bio_v1) .+ offsets["si1"]))
        push!(ymins, minimum(Float64.(bio_v2) .+ offsets["si2"]))
        push!(ymaxs, maximum(Float64.(bio_v2) .+ offsets["si2"]))
    end
    for mode_name in MODE_LABELS
        tr = traces[mode_name]
        isempty(tr) && continue
        append!(xs, [first(tr.time_s), last(tr.time_s)])
        vals = Float64.(tr.V1) .+ offsets[mode_name]
        push!(ymins, minimum(vals))
        push!(ymaxs, maximum(vals))
        if show_si3 && :V3 in propertynames(tr)
            si3_vals = Float64.(tr.V3)
            finite = isfinite.(si3_vals)
            if any(finite)
                vals = si3_vals[finite] .+ offsets["$(mode_name)_si3"]
                push!(ymins, minimum(vals))
                push!(ymaxs, maximum(vals))
            end
        end
    end
    sm_trace = get(traces, "presynaptic", DataFrame())
    if isempty(sm_trace) || !(:sm in propertynames(sm_trace)) || !any(isfinite, sm_trace.sm)
        sm_trace = get(traces, "postsynaptic", DataFrame())
    end
    if !isempty(sm_trace) && :sm in propertynames(sm_trace)
        _, sm_vals = neuromod_trace_vectors(sm_trace.time_s, sm_trace.sm, shared_neuromod_offset(offsets))
        if !isempty(sm_vals)
            push!(ymins, minimum(sm_vals))
            push!(ymaxs, maximum(sm_vals))
        end
    end
    return minimum(xs), maximum(xs), minimum(ymins) - 12.0, maximum(ymaxs) + 60.0
end

function reset_trace_limits!(ax::Axis, bio_t1, bio_v1, bio_t2, bio_v2, traces::Dict{String, DataFrame}; source::Symbol = :representative, show_si3::Bool = false)
    x_min, x_max, y_min, y_max = trace_panel_limits(bio_t1, bio_v1, bio_t2, bio_v2, traces; source = source, show_si3 = show_si3)
    xlims!(ax, x_min, x_max)
    ylims!(ax, y_min, y_max)
end

function set_model_trace_observables!(trace_obs, traces::Dict{String, DataFrame})
    for mode_name in MODE_LABELS
        tr = traces[mode_name]
        if isempty(tr)
            trace_obs[(mode_name, :time)][] = Float64[]
            trace_obs[(mode_name, :V0)][] = Float64[]
            trace_obs[(mode_name, :V1)][] = Float64[]
            trace_obs[(mode_name, :V3)][] = Float64[]
            trace_obs[(mode_name, :sm)][] = Float64[]
        else
            trace_obs[(mode_name, :time)][] = Float64.(tr.time_s)
            trace_obs[(mode_name, :V0)][] = Float64.(tr.V0)
            trace_obs[(mode_name, :V1)][] = Float64.(tr.V1)
            trace_obs[(mode_name, :V3)][] = :V3 in propertynames(tr) ? Float64.(tr.V3) : fill(NaN, nrow(tr))
            trace_obs[(mode_name, :sm)][] = :sm in propertynames(tr) ? Float64.(tr.sm) : fill(NaN, nrow(tr))
        end
    end
end

function traces_have_finite_si3(traces::Dict{String, DataFrame})
    for mode_name in MODE_LABELS
        tr = traces[mode_name]
        if isempty(tr) || !(:V3 in propertynames(tr)) || !any(isfinite, tr.V3)
            return false
        end
    end
    return true
end

function selected_gain_vectors(summary::DataFrame, selected::Tuple{String, Float64}, ycol::Symbol)
    row = row_for_selection(summary, selected[1], selected[2])
    row === nothing && return Float64[], Float64[]
    y = row[ycol]
    (!isfinite(Float64(y))) && return Float64[], Float64[]
    return [Float64(row.control_gain)], [Float64(y)]
end

function matched_trace_vectors(time_s::AbstractVector{<:Real}, voltage_mv::AbstractVector{<:Real}, offset::Float64)
    n = min(length(time_s), length(voltage_mv))
    n == 0 && return Float64[], Float64[]
    return Float64.(time_s[1:n]), Float64.(voltage_mv[1:n]) .+ offset
end

safe_x(x_obs, y_obs) = lift((x, y) -> matched_xy_vectors(x, y)[1], x_obs, y_obs)
safe_y(x_obs, y_obs) = lift((x, y) -> matched_xy_vectors(x, y)[2], x_obs, y_obs)
trend_x_obs(x_obs, y_obs) = lift((x, y) -> trend_line_vectors(matched_xy_vectors(x, y)...)[1], x_obs, y_obs)
trend_y_obs(x_obs, y_obs) = lift((x, y) -> trend_line_vectors(matched_xy_vectors(x, y)...)[2], x_obs, y_obs)

function add_labeled_checkbox!(grid, row::Int, label::String; checked::Bool = true)
    cb = Checkbox(grid[row, 1]; checked = checked)
    Label(grid[row, 2], label; halign = :left, tellwidth = false)
    return cb
end

function refresh_scan_observables!(summary_obs, raw_obs, plot_obs)
    summary = summary_obs[]
    raw = raw_obs[]
    for mode_name in MODE_LABELS
        on_x, on_y = summary_vectors(summary, mode_name, :mean_pre_frequency_hz)
        off_x, off_y = summary_vectors(summary, mode_name, :mean_post_frequency_hz)
        ron_x, ron_y = raw_vectors(raw, mode_name, "pre")
        roff_x, roff_y = raw_vectors(raw, mode_name, "post")
        plot_obs[(mode_name, :on_summary_x)][] = on_x
        plot_obs[(mode_name, :on_summary_y)][] = on_y
        plot_obs[(mode_name, :off_summary_x)][] = off_x
        plot_obs[(mode_name, :off_summary_y)][] = off_y
        plot_obs[(mode_name, :on_raw_x)][] = ron_x
        plot_obs[(mode_name, :on_raw_y)][] = ron_y
        plot_obs[(mode_name, :off_raw_x)][] = roff_x
        plot_obs[(mode_name, :off_raw_y)][] = roff_y
    end
end

function plot_gain_panel!(ax::Axis, summary::DataFrame, raw::DataFrame, phase::String)
    ycol = phase == "pre" ? :mean_pre_frequency_hz : :mean_post_frequency_hz
    for mode_name in MODE_LABELS
        color = DASH_COLORS[mode_name]
        marker = DASH_MARKERS[mode_name]
        rx, ry = raw_vectors(raw, mode_name, phase)
        !isempty(rx) && scatter!(ax, rx, ry; color = (color, 0.22), markersize = 6, strokewidth = 0)
        x, y = summary_vectors(summary, mode_name, ycol)
        scatter!(ax, x, y; color = color, marker = marker, markersize = 10, strokecolor = :black, strokewidth = 0.8)
    end
    ax.xlabel = "Gain"
    ax.ylabel = phase == "pre" ? "Freq on (Hz)" : "Freq off (Hz)"
    ax.title = phase == "pre" ? "During Drive" : "After Drive Removal"
    xlims!(ax, 0, 8)
    ylims!(ax, phase == "post" ? -0.005 : 0.0, 0.22)
end

function plot_representative_traces!(ax::Axis, bio_t1, bio_v1, bio_t2, bio_v2, traces::Dict{String, DataFrame}; source::Symbol = :representative, show_si3::Bool = false, show_burst_markers::Bool = true, detector = (spike_threshold_mv = -20.0, spike_refractory_s = 0.02, burst_factor = 3.0))
    offsets = trace_offsets()
    if source == :scan_protocol && !isempty(traces["presynaptic"]) && !isempty(traces["postsynaptic"])
        lines!(ax, traces["presynaptic"].time_s, Float64.(traces["presynaptic"].V0) .+ offsets["si1"]; color = RGBf(0.35, 0.35, 0.35), linewidth = 1.4)
    else
        lines!(ax, bio_t1, Float64.(bio_v1) .+ offsets["si1"]; color = RGBf(0.35, 0.35, 0.35), linewidth = 1.4)
        lines!(ax, bio_t2, Float64.(bio_v2) .+ offsets["si2"]; color = RGBf(0.1, 0.1, 0.1), linewidth = 1.4)
        show_burst_markers && plot_burst_markers!(ax, bio_t2, bio_v2, offsets["si2"]; color = :black, marker = :utriangle, strokecolor = :white, detector = detector)
    end
    for mode_name in MODE_LABELS
        tr = traces[mode_name]
        isempty(tr) && continue
        lines!(ax, tr.time_s, Float64.(tr.V1) .+ offsets[mode_name]; color = DASH_COLORS[mode_name], linewidth = 1.4)
        show_burst_markers && plot_burst_markers!(ax, tr.time_s, tr.V1, offsets[mode_name]; color = DASH_COLORS[mode_name], marker = DASH_MARKERS[mode_name], detector = detector)
        if show_si3 && :V3 in propertynames(tr)
            lines!(ax, tr.time_s, Float64.(tr.V3) .+ offsets["$(mode_name)_si3"]; color = DASH_COLORS[mode_name], linewidth = 1.4)
            show_burst_markers && plot_burst_markers!(ax, tr.time_s, tr.V3, offsets["$(mode_name)_si3"]; color = DASH_COLORS[mode_name], marker = DASH_MARKERS[mode_name], detector = detector)
        end
    end
    sm_trace = !isempty(traces["presynaptic"]) ? traces["presynaptic"] : traces["postsynaptic"]
    if !isempty(sm_trace) && :sm in propertynames(sm_trace)
        sm_x, sm_y = neuromod_trace_vectors(sm_trace.time_s, sm_trace.sm, shared_neuromod_offset(offsets))
        lines!(ax, sm_x, sm_y; color = RGBf(0.45, 0.24, 0.62), linewidth = 1.1)
    end
    x_min, x_max, y_min, y_max = trace_panel_limits(bio_t1, bio_v1, bio_t2, bio_v2, traces; source = source, show_si3 = show_si3)
    label_x = x_min + 0.045 * (x_max - x_min)
    if source == :scan_protocol
        text!(ax, label_x, offsets["si1"] + 38; text = "Si1 presyn scan", align = (:left, :center), fontsize = 15, color = RGBf(0.35, 0.35, 0.35))
    else
        text!(ax, label_x, offsets["si1"] + 38; text = "Si1 biology", align = (:left, :center), fontsize = 15, color = RGBf(0.35, 0.35, 0.35))
        text!(ax, label_x, offsets["si2"] + 38; text = "Si2 biology", align = (:left, :center), fontsize = 15, color = RGBf(0.1, 0.1, 0.1))
    end
    text!(ax, label_x, offsets["presynaptic"] + 38; text = "Si2 presynaptic", align = (:left, :center), fontsize = 15, color = DASH_COLORS["presynaptic"])
    if show_si3
        text!(ax, label_x, offsets["presynaptic_si3"] + 38; text = "Si3 presynaptic", align = (:left, :center), fontsize = 15, color = DASH_COLORS["presynaptic"])
    end
    text!(ax, label_x, offsets["postsynaptic"] + 38; text = "Si2 postsynaptic", align = (:left, :center), fontsize = 15, color = DASH_COLORS["postsynaptic"])
    text!(ax, label_x, shared_neuromod_offset(offsets) + 22; text = "sm neuromodulation", align = (:left, :center), fontsize = 12, color = RGBf(0.45, 0.24, 0.62))
    if show_si3
        text!(ax, label_x, offsets["postsynaptic_si3"] + 38; text = "Si3 postsynaptic", align = (:left, :center), fontsize = 15, color = DASH_COLORS["postsynaptic"])
    end
    ax.title = "Representative Traces"
    ax.ylabel = "V (mV)"
    ax.xticksvisible = false
    ax.yticksvisible = false
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    hidespines!(ax, :t, :r)
    xlims!(ax, x_min, x_max)
    ylims!(ax, y_min, y_max)
    add_l_scale_bar!(ax, x_min, x_max, y_min, y_max)
end

function build_clean_figure(summary::DataFrame, raw::DataFrame, traces::Dict{String, DataFrame}; trace_source::Symbol = :representative, show_si3::Bool = false, show_burst_markers::Bool = true, detector = (spike_threshold_mv = -20.0, spike_refractory_s = 0.02, burst_factor = 3.0), time_scale::Float64 = REPRESENTATIVE_BIOLOGY_TIME_SCALE)
    bio_t1, bio_v1, bio_t2, bio_v2 = load_biology_pair(; time_scale = time_scale)
    fig = Figure(size = (1180, 1320), backgroundcolor = :white, fontsize = 18)
    axA = Axis(fig[1, 1], backgroundcolor = RGBf(0.98, 0.98, 0.99))
    axB = Axis(fig[1, 2], backgroundcolor = RGBf(0.98, 0.98, 0.99))
    axC = Axis(fig[2, 1], backgroundcolor = RGBf(0.98, 0.98, 0.99))
    axD = Axis(fig[2, 2], backgroundcolor = RGBf(0.98, 0.98, 0.99))
    axE = Axis(fig[3, 1:2], backgroundcolor = RGBf(0.98, 0.98, 0.99))
    plot_circuit_panel!(axA)
    plot_compare_panel!(axB, representative_model_points(traces; detector_kwargs(detector)...); detector = detector, time_scale = time_scale)
    plot_gain_panel!(axC, summary, raw, "pre")
    plot_gain_panel!(axD, summary, raw, "post")
    plot_representative_traces!(axE, bio_t1, bio_v1, bio_t2, bio_v2, traces; source = trace_source, show_si3 = show_si3, show_burst_markers = show_burst_markers, detector = detector)
    add_panel_label!(axA, "A")
    add_panel_label!(axB, "B")
    add_panel_label!(axC, "C")
    add_panel_label!(axD, "D")
    add_panel_label!(axE, "E")
    rowsize!(fig.layout, 1, Fixed(235))
    rowsize!(fig.layout, 2, Fixed(335))
    rowsize!(fig.layout, 3, Fixed(255))
    rowgap!(fig.layout, 18)
    colgap!(fig.layout, 16)
    return fig
end

function build_dashboard()
    raw_obs = Observable(load_dashboard_raw())
    summary_obs = Observable(load_dashboard_summary())
    selected_obs = Observable(("presynaptic", PRESYNAPTIC_ANCHOR_GAIN))
    traces_obs = Observable(load_representative_traces())
    full_initial_state = Ref{Union{Nothing, Vector{Float64}}}(nothing)
    detector_obs = Observable((spike_threshold_mv = -20.0, spike_refractory_s = 0.02, burst_factor = 3.0))
    trace_time_scale_obs = Observable(REPRESENTATIVE_BIOLOGY_TIME_SCALE)
    model_points_obs = Observable(representative_model_points(traces_obs[]; detector_kwargs(detector_obs[])...))
    bio_points_obs = Observable(load_biological_trace_points(; time_scale = trace_time_scale_obs[], detector_kwargs(detector_obs[])...))
    status_obs = Observable("Ready")

    plot_obs = Dict{Tuple{String, Symbol}, Observable{Vector{Float64}}}()
    for mode_name in MODE_LABELS, key in [:on_summary_x, :on_summary_y, :off_summary_x, :off_summary_y, :on_raw_x, :on_raw_y, :off_raw_x, :off_raw_y]
        plot_obs[(mode_name, key)] = Observable(Float64[])
    end
    refresh_scan_observables!(summary_obs, raw_obs, plot_obs)

    point_obs = Dict{Tuple{String, Symbol}, Observable{Vector{Float64}}}()
    for mode_name in MODE_LABELS, key in [:x, :y]
        point_obs[(mode_name, key)] = Observable(Float64[])
    end
    refresh_model_point_observables!(model_points_obs, point_obs)

    trace_obs = Dict{Tuple{String, Symbol}, Observable{Vector{Float64}}}()
    for mode_name in MODE_LABELS
        trace_obs[(mode_name, :time)] = Observable(Float64[])
        trace_obs[(mode_name, :V0)] = Observable(Float64[])
        trace_obs[(mode_name, :V1)] = Observable(Float64[])
        trace_obs[(mode_name, :V3)] = Observable(Float64[])
        trace_obs[(mode_name, :sm)] = Observable(Float64[])
    end
    set_model_trace_observables!(trace_obs, traces_obs[])
    bio_t1, bio_v1, bio_t2, bio_v2 = load_biology_pair(; time_scale = trace_time_scale_obs[])
    bio_t1_obs = Observable(bio_t1)
    bio_v1_obs = Observable(bio_v1)
    bio_t2_obs = Observable(bio_t2)
    bio_v2_obs = Observable(bio_v2)
    lit_x, lit_y = load_biological_fig2c(; time_scale = trace_time_scale_obs[])
    lit_x_obs = Observable(lit_x)
    lit_y_obs = Observable(lit_y)
    fig = Figure(size = (1980, 1380), backgroundcolor = :white, fontsize = 15)
    controls = GridLayout(fig[1:2, 1], tellwidth = true)
    parameter_controls = GridLayout(fig[1:2, 2], tellwidth = true)
    plots = GridLayout(fig[1:2, 3:4])

    Label(controls[1, 1:2], "Scan Controls"; fontsize = 20, font = :bold, halign = :left)
    cb_presyn_scan = add_labeled_checkbox!(controls, 2, "scan presynaptic"; checked = true)
    cb_postsyn_scan = add_labeled_checkbox!(controls, 3, "scan postsynaptic"; checked = true)
    order_menu = Menu(controls[4, 1:2]; options = ORDER_LABELS, default = "left-to-right", width = 190)
    Label(controls[5, 1:2], "scan order"; halign = :left)

    sg_scan = SliderGrid(
        controls[6, 1:2],
        (label = "min gain", range = 0.0:0.1:8.0, startvalue = 0.0, format = "{:.1f}"),
        (label = "max gain", range = 0.5:0.1:10.0, startvalue = 8.0, format = "{:.1f}"),
        (label = "points", range = 3:1:120, startvalue = 36, format = x -> string(Int(round(x)))),
        (label = "max time s", range = 120.0:60.0:2400.0, startvalue = 900.0, format = "{:.0f}"),
        (label = "saveat ms", range = 2.0:1.0:25.0, startvalue = 5.0, format = "{:.0f}"),
    )
    run_scan_button = Button(controls[7, 1:2]; label = "Run Scan", width = 140)

    Label(controls[8, 1:2], "Trace Controls"; fontsize = 20, font = :bold, halign = :left)
    cb_trace_scan_protocol = add_labeled_checkbox!(controls, 9, "use scan burst protocol"; checked = false)
    cb_show_si3 = add_labeled_checkbox!(controls, 10, "show Si3 traces"; checked = true)
    cb_show_burst_markers = add_labeled_checkbox!(controls, 11, "show burst onsets"; checked = true)
    sg_timescale = SliderGrid(
        controls[12, 1:2],
        (label = "time scale", range = 1.0:0.25:10.0, startvalue = trace_time_scale_obs[], format = x -> string(round(Float64(x), digits = 2), "x")),
    )
    sg_trace = SliderGrid(
        controls[13, 1:2],
        (label = "presyn gain", range = 0.0:0.05:10.0, startvalue = PRESYNAPTIC_ANCHOR_GAIN, format = "{:.2f}"),
        (label = "postsyn gain", range = 0.0:0.05:10.0, startvalue = POSTSYNAPTIC_ANCHOR_GAIN, format = "{:.2f}"),
        (label = "Ca0 IC", range = 0.0:0.05:2.0, startvalue = initial_state()[11], format = "{:.2f}"),
        (label = "Ca1 IC", range = 0.0:0.05:2.0, startvalue = initial_state()[12], format = "{:.2f}"),
        (label = "Ca2 IC", range = 0.0:0.05:2.0, startvalue = initial_state()[13], format = "{:.2f}"),
        (label = "Ca3 IC", range = 0.0:0.05:2.0, startvalue = initial_state()[14], format = "{:.2f}"),
        (label = "Ca4 IC", range = 0.0:0.05:2.0, startvalue = initial_state()[15], format = "{:.2f}"),
        (label = "spike threshold", range = -50.0:1.0:0.0, startvalue = detector_obs[].spike_threshold_mv, format = "{:.0f}"),
        (label = "spike refractory s", range = 0.005:0.005:0.100, startvalue = detector_obs[].spike_refractory_s, format = "{:.3f}"),
        (label = "burst factor", range = 1.1:0.1:5.0, startvalue = detector_obs[].burst_factor, format = "{:.1f}"),
    )

    terminal_ic_button = Button(controls[14, 1:2]; label = "Set ICs to Terminal State", width = 220)
    save_button = Button(controls[15, 1]; label = "Save Clean PNG", width = 140)
    write_csv_button = Button(controls[15, 2]; label = "Write CSVs", width = 120)
    Label(controls[16, 1:2], status_obs; halign = :left, tellwidth = false)

    Label(parameter_controls[1, 1:2], "Parameters"; fontsize = 20, font = :bold, halign = :left)
    param_defaults = load_param_defaults()
    param_defaults_by_key = Dict(key => param_defaults[i] for (i, key) in enumerate(DASHBOARD_PARAM_KEYS))
    param_specs = Dict{Symbol, Any}(
        :x_shift_si2 => (label = "x_shift1/2", range = -8.0:0.05:-4.0, startvalue = param_defaults_by_key[:x_shift_si2], format = "{:.2f}"),
        :x_shift_si3 => (label = "x_shift3/4", range = -8.0:0.05:-4.0, startvalue = param_defaults_by_key[:x_shift_si3], format = "{:.2f}"),
        :presyn_g0 => (label = "g0 presyn", range = 0.0:0.0001:0.0100, startvalue = param_defaults_by_key[:presyn_g0], format = "{:.4f}"),
        :postsyn_g0 => (label = "g0 postsyn", range = 0.0:0.0001:0.0100, startvalue = param_defaults_by_key[:postsyn_g0], format = "{:.4f}"),
        :alpha1 => (label = "alpha1", range = 0.001:0.001:0.050, startvalue = param_defaults_by_key[:alpha1], format = "{:.3f}"),
        :beta1 => (label = "beta1", range = 0.0001:0.0001:0.0100, startvalue = param_defaults_by_key[:beta1], format = "{:.4f}"),
        :s0_floor => (label = "s0 floor", range = 0.01:0.0025:0.5000, startvalue = param_defaults_by_key[:s0_floor], format = "{:.4f}"),
        :alpham => (label = "alpham", range = 0.0005:0.0005:0.1000, startvalue = param_defaults_by_key[:alpham], format = "{:.4f}"),
        :betam => (label = "betam", range = 0.00005:0.00005:0.00200, startvalue = param_defaults_by_key[:betam], format = "{:.5f}"),
        :sm_floor => (label = "sm floor", range = 0.01:0.0025:0.5000, startvalue = param_defaults_by_key[:sm_floor], format = "{:.4f}"),
        :ca_shift_si2 => (label = "Ca_shift1/2", range = -120.0:1.0:20.0, startvalue = param_defaults_by_key[:ca_shift_si2], format = "{:.0f}"),
        :ca_shift_si3 => (label = "Ca_shift3/4", range = -120.0:1.0:20.0, startvalue = param_defaults_by_key[:ca_shift_si3], format = "{:.0f}"),
        :presynaptic_base_g => (label = "presyn g41/g32 base", range = 0.0:0.0001:0.0100, startvalue = param_defaults_by_key[:presynaptic_base_g], format = "{:.4f}"),
        :direct_post_base_g => (label = "postsyn g41/g32 base", range = 0.0:0.0001:0.0100, startvalue = param_defaults_by_key[:direct_post_base_g], format = "{:.4f}"),
        :presyn_alphax => (label = "presyn alphax", range = 0.001:0.001:0.050, startvalue = param_defaults_by_key[:presyn_alphax], format = "{:.3f}"),
        :postsyn_alphax => (label = "postsyn alphax", range = 0.001:0.001:0.100, startvalue = param_defaults_by_key[:postsyn_alphax], format = "{:.3f}"),
        :presyn_betax => (label = "presyn betax", range = 0.0001:0.0001:0.0100, startvalue = param_defaults_by_key[:presyn_betax], format = "{:.4f}"),
        :postsyn_betax => (label = "postsyn betax", range = 0.0001:0.0001:0.0100, startvalue = param_defaults_by_key[:postsyn_betax], format = "{:.4f}"),
        :si3_exc_floor => (label = "s3/s4 floor", range = 0.01:0.0025:0.5000, startvalue = param_defaults_by_key[:si3_exc_floor], format = "{:.4f}"),
        :slow_inhib_g => (label = "g14/g23 slow inhib", range = 0.0:0.0001:0.0300, startvalue = param_defaults_by_key[:slow_inhib_g], format = "{:.4f}"),
        :slow_inhib_alpha => (label = "alphai", range = 0.001:0.001:0.050, startvalue = param_defaults_by_key[:slow_inhib_alpha], format = "{:.3f}"),
        :slow_inhib_beta => (label = "betai", range = 0.0001:0.0001:0.0200, startvalue = param_defaults_by_key[:slow_inhib_beta], format = "{:.4f}"),
        :si2_mutual_inhib_g => (label = "Si2 mutual inhib g", range = 0.0:0.0001:0.0300, startvalue = param_defaults_by_key[:si2_mutual_inhib_g], format = "{:.4f}"),
        :si2_mutual_inhib_alpha => (label = "Si2 mutual inhib alpha", range = 0.001:0.001:0.050, startvalue = param_defaults_by_key[:si2_mutual_inhib_alpha], format = "{:.3f}"),
        :si2_mutual_inhib_beta => (label = "Si2 mutual inhib beta", range = 0.0001:0.0001:0.0200, startvalue = param_defaults_by_key[:si2_mutual_inhib_beta], format = "{:.4f}"),
        :si_inhib_floor => (label = "s1/s2 floor", range = 0.01:0.0025:0.5000, startvalue = param_defaults_by_key[:si_inhib_floor], format = "{:.4f}"),
        :si3_mutual_inhib_g => (label = "Si3 mutual inhib g", range = 0.0:0.0001:0.0300, startvalue = param_defaults_by_key[:si3_mutual_inhib_g], format = "{:.4f}"),
        :si3_mutual_inhib_alpha => (label = "Si3 mutual inhib alpha", range = 0.001:0.001:0.050, startvalue = param_defaults_by_key[:si3_mutual_inhib_alpha], format = "{:.3f}"),
        :si3_mutual_inhib_beta => (label = "Si3 mutual inhib beta", range = 0.0001:0.0001:0.0200, startvalue = param_defaults_by_key[:si3_mutual_inhib_beta], format = "{:.4f}"),
        :t1_ms => (label = "t1", range = 1000.0:1000.0:100000.0, startvalue = param_defaults_by_key[:t1_ms], format = "{:.0f}"),
    )
    param_sections = [
        ("Voltage/Ca Shifts", [:x_shift_si2, :x_shift_si3, :ca_shift_si2, :ca_shift_si3]),
        ("Si1 Excitation", [:presyn_g0, :postsyn_g0, :alpha1, :beta1, :s0_floor]),
        ("Modulation/Timing", [:alpham, :betam, :sm_floor, :t1_ms]),
        ("Excitatory Coupling", [:presynaptic_base_g, :direct_post_base_g, :presyn_alphax, :postsyn_alphax, :presyn_betax, :postsyn_betax, :si3_exc_floor]),
        ("Slow Inhibition", [:slow_inhib_g, :slow_inhib_alpha, :slow_inhib_beta, :si_inhib_floor]),
        ("Mutual Inhibition", [:si2_mutual_inhib_g, :si2_mutual_inhib_alpha, :si2_mutual_inhib_beta, :si3_mutual_inhib_g, :si3_mutual_inhib_alpha, :si3_mutual_inhib_beta]),
    ]
    param_grids = Any[]
    param_slider_by_key = Dict{Symbol, Any}()
    row = 2
    for (title, keys) in param_sections
        Label(parameter_controls[row, 1:2], title; fontsize = 16, font = :bold, halign = :left)
        grid = SliderGrid(parameter_controls[row + 1, 1:2], (param_specs[key] for key in keys)...)
        push!(param_grids, grid)
        for (key, slider) in zip(keys, grid.sliders)
            param_slider_by_key[key] = slider
        end
        row += 2
    end
    sg_params = (sliders = reduce(vcat, [collect(grid.sliders) for grid in param_grids]), grids = param_grids, by_key = param_slider_by_key)

    active_param_text = Observable("Active parameters pending")
    Label(parameter_controls[row, 1:2], active_param_text; fontsize = 11, halign = :left, tellwidth = false)
    row += 1

    match_baselines_button = Button(parameter_controls[row, 1:2]; label = "Match Post Baseline to Pre", width = 290)
    row += 1

    reset_params_button = Button(parameter_controls[row, 1]; label = "Reset Params", width = 140)
    save_params_button = Button(parameter_controls[row, 2]; label = "Save Params", width = 140)

    axA = Axis(plots[1, 1], backgroundcolor = RGBf(0.98, 0.98, 0.99))
    plot_circuit_panel!(axA)
    axB = Axis(plots[1, 2], title = "Si1 Rate vs Burst Rate", xlabel = "Si1 spike freq (Hz)", ylabel = "Burst freq (Hz)", backgroundcolor = RGBf(0.98, 0.98, 0.99))
    bio_x = lift(points -> finite_cols(points, :si1_frequency_hz, :burst_frequency_hz)[1], bio_points_obs)
    bio_y = lift(points -> finite_cols(points, :si1_frequency_hz, :burst_frequency_hz)[2], bio_points_obs)
    scatter!(axB, safe_x(lit_x_obs, lit_y_obs), safe_y(lit_x_obs, lit_y_obs); color = LITERATURE_GREY, marker = :circle, markersize = 10, strokecolor = :black, strokewidth = 0.8)
    scatter!(axB, safe_x(bio_x, bio_y), safe_y(bio_x, bio_y); color = :black, marker = :utriangle, markersize = 10, strokecolor = :white, strokewidth = 0.8)
    trend_y_observables = Observable{Vector{Float64}}[]
    for mode_name in MODE_LABELS
        trend_y = trend_y_obs(point_obs[(mode_name, :x)], point_obs[(mode_name, :y)])
        push!(trend_y_observables, trend_y)
        scatter!(axB, safe_x(point_obs[(mode_name, :x)], point_obs[(mode_name, :y)]), safe_y(point_obs[(mode_name, :x)], point_obs[(mode_name, :y)]);
            color = (DASH_COLORS[mode_name], 0.84), marker = DASH_MARKERS[mode_name], markersize = 10, strokecolor = :black, strokewidth = 0.8)
    end
    lines!(axB, trend_x_obs(lit_x_obs, lit_y_obs), trend_y_obs(lit_x_obs, lit_y_obs); color = LITERATURE_GREY, linewidth = 3.0)
    bio_trend_x = trend_x_obs(bio_x, bio_y)
    bio_trend_y = trend_y_obs(bio_x, bio_y)
    lines!(axB, bio_trend_x, bio_trend_y; color = :black, linewidth = 4.0)
    for (i, mode_name) in enumerate(MODE_LABELS)
        trend_x = trend_x_obs(point_obs[(mode_name, :x)], point_obs[(mode_name, :y)])
        trend_y = trend_y_observables[i]
        lines!(axB, trend_x, trend_y; color = DASH_COLORS[mode_name], linewidth = 4.0)
    end
    xlims!(axB, 0, 6.5)
    function refresh_compare_ylim!()
        ylims!(axB, 0, compare_panel_ylimit(model_points_obs[], bio_points_obs[]))
    end
    on(model_points_obs) do _
        refresh_compare_ylim!()
    end
    on(bio_points_obs) do _
        refresh_compare_ylim!()
    end
    refresh_compare_ylim!()
    hidespines!(axB, :t, :r)
    add_panel_label!(axA, "A")
    add_panel_label!(axB, "B")

    axC = Axis(plots[2, 1], title = "During Drive", xlabel = "Gain", ylabel = "Freq on (Hz)", backgroundcolor = RGBf(0.98, 0.98, 0.99))
    axD = Axis(plots[2, 2], title = "After Drive Removal", xlabel = "Gain", ylabel = "Freq off (Hz)", backgroundcolor = RGBf(0.98, 0.98, 0.99))
    pick_targets = IdDict{Any, Tuple{String, Symbol}}()
    for mode_name in MODE_LABELS
        color = DASH_COLORS[mode_name]
        marker = DASH_MARKERS[mode_name]
        scatter!(axC, safe_x(plot_obs[(mode_name, :on_raw_x)], plot_obs[(mode_name, :on_raw_y)]), safe_y(plot_obs[(mode_name, :on_raw_x)], plot_obs[(mode_name, :on_raw_y)]); color = (color, 0.20), markersize = 6, strokewidth = 0)
        sc_on = scatter!(axC, safe_x(plot_obs[(mode_name, :on_summary_x)], plot_obs[(mode_name, :on_summary_y)]), safe_y(plot_obs[(mode_name, :on_summary_x)], plot_obs[(mode_name, :on_summary_y)]); color = color, marker = marker, markersize = 11, strokecolor = :black, strokewidth = 0.8)
        pick_targets[sc_on] = (mode_name, :mean_pre_frequency_hz)
        scatter!(axD, safe_x(plot_obs[(mode_name, :off_raw_x)], plot_obs[(mode_name, :off_raw_y)]), safe_y(plot_obs[(mode_name, :off_raw_x)], plot_obs[(mode_name, :off_raw_y)]); color = (color, 0.20), markersize = 6, strokewidth = 0)
        sc_off = scatter!(axD, safe_x(plot_obs[(mode_name, :off_summary_x)], plot_obs[(mode_name, :off_summary_y)]), safe_y(plot_obs[(mode_name, :off_summary_x)], plot_obs[(mode_name, :off_summary_y)]); color = color, marker = marker, markersize = 11, strokecolor = :black, strokewidth = 0.8)
        pick_targets[sc_off] = (mode_name, :mean_post_frequency_hz)
    end
    selected_on_x = lift((summary, selected) -> selected_gain_vectors(summary, selected, :mean_pre_frequency_hz)[1], summary_obs, selected_obs)
    selected_on_y = lift((summary, selected) -> selected_gain_vectors(summary, selected, :mean_pre_frequency_hz)[2], summary_obs, selected_obs)
    selected_off_x = lift((summary, selected) -> selected_gain_vectors(summary, selected, :mean_post_frequency_hz)[1], summary_obs, selected_obs)
    selected_off_y = lift((summary, selected) -> selected_gain_vectors(summary, selected, :mean_post_frequency_hz)[2], summary_obs, selected_obs)
    scatter!(axC, safe_x(selected_on_x, selected_on_y), safe_y(selected_on_x, selected_on_y); color = RGBf(1.0, 0.78, 0.08), marker = :star5, markersize = 22, strokecolor = :black, strokewidth = 1.0)
    scatter!(axD, safe_x(selected_off_x, selected_off_y), safe_y(selected_off_x, selected_off_y); color = RGBf(1.0, 0.78, 0.08), marker = :star5, markersize = 22, strokecolor = :black, strokewidth = 1.0)
    for ax in (axC, axD)
        xlims!(ax, 0, 8)
        ylims!(ax, ax === axD ? -0.005 : 0.0, 0.22)
    end
    add_panel_label!(axC, "C")
    add_panel_label!(axD, "D")

    axE = Axis(plots[3, 1:2], backgroundcolor = RGBf(0.98, 0.98, 0.99))
    offsets = trace_offsets()
    scan_protocol_visible = cb_trace_scan_protocol.checked
    empirical_protocol_visible = lift(x -> !x, scan_protocol_visible)
    empirical_burst_markers_visible = lift((empirical, show_markers) -> empirical && show_markers, empirical_protocol_visible, cb_show_burst_markers.checked)
    bio_si1_y = lift(v -> Float64.(v) .+ offsets["si1"], bio_v1_obs)
    bio_si2_y = lift(v -> Float64.(v) .+ offsets["si2"], bio_v2_obs)
    lines!(axE, bio_t1_obs, bio_si1_y; color = RGBf(0.35, 0.35, 0.35), linewidth = 1.4, visible = empirical_protocol_visible)
    lines!(axE, bio_t2_obs, bio_si2_y; color = RGBf(0.1, 0.1, 0.1), linewidth = 1.4, visible = empirical_protocol_visible)
    bio_tick_x = lift((t, v, d) -> burst_tick_vectors(t, v, offsets["si2"]; detector_kwargs(d)...)[1], bio_t2_obs, bio_v2_obs, detector_obs)
    bio_tick_y = lift((t, v, d) -> burst_tick_vectors(t, v, offsets["si2"]; detector_kwargs(d)...)[2], bio_t2_obs, bio_v2_obs, detector_obs)
    bio_marker_x = lift((t, v, d) -> burst_marker_vectors(t, v, offsets["si2"]; detector_kwargs(d)...)[1], bio_t2_obs, bio_v2_obs, detector_obs)
    bio_marker_y = lift((t, v, d) -> burst_marker_vectors(t, v, offsets["si2"]; detector_kwargs(d)...)[2], bio_t2_obs, bio_v2_obs, detector_obs)
    lines!(axE, bio_tick_x, bio_tick_y; color = (:black, 0.28), linewidth = 1.0, visible = empirical_burst_markers_visible)
    scatter!(axE, bio_marker_x, bio_marker_y; color = :black, marker = :utriangle, markersize = 10, strokecolor = :white, strokewidth = 0.8, visible = empirical_burst_markers_visible)
    pre_drive_x = lift((t, v) -> matched_trace_vectors(t, v, offsets["si1"])[1], trace_obs[("presynaptic", :time)], trace_obs[("presynaptic", :V0)])
    pre_drive_y = lift((t, v) -> matched_trace_vectors(t, v, offsets["si1"])[2], trace_obs[("presynaptic", :time)], trace_obs[("presynaptic", :V0)])
    pre_si2_x = lift((t, v) -> matched_trace_vectors(t, v, offsets["presynaptic"])[1], trace_obs[("presynaptic", :time)], trace_obs[("presynaptic", :V1)])
    pre_si2_y = lift((t, v) -> matched_trace_vectors(t, v, offsets["presynaptic"])[2], trace_obs[("presynaptic", :time)], trace_obs[("presynaptic", :V1)])
    post_si2_x = lift((t, v) -> matched_trace_vectors(t, v, offsets["postsynaptic"])[1], trace_obs[("postsynaptic", :time)], trace_obs[("postsynaptic", :V1)])
    post_si2_y = lift((t, v) -> matched_trace_vectors(t, v, offsets["postsynaptic"])[2], trace_obs[("postsynaptic", :time)], trace_obs[("postsynaptic", :V1)])
    shared_sm_x = lift((pre_t, pre_sm, post_t, post_sm) -> begin
        any(isfinite, pre_sm) && return neuromod_trace_vectors(pre_t, pre_sm, shared_neuromod_offset(offsets))[1]
        neuromod_trace_vectors(post_t, post_sm, shared_neuromod_offset(offsets))[1]
    end, trace_obs[("presynaptic", :time)], trace_obs[("presynaptic", :sm)], trace_obs[("postsynaptic", :time)], trace_obs[("postsynaptic", :sm)])
    shared_sm_y = lift((pre_t, pre_sm, post_t, post_sm) -> begin
        any(isfinite, pre_sm) && return neuromod_trace_vectors(pre_t, pre_sm, shared_neuromod_offset(offsets))[2]
        neuromod_trace_vectors(post_t, post_sm, shared_neuromod_offset(offsets))[2]
    end, trace_obs[("presynaptic", :time)], trace_obs[("presynaptic", :sm)], trace_obs[("postsynaptic", :time)], trace_obs[("postsynaptic", :sm)])
    pre_si3_x = lift((t, v) -> matched_trace_vectors(t, v, offsets["presynaptic_si3"])[1], trace_obs[("presynaptic", :time)], trace_obs[("presynaptic", :V3)])
    pre_si3_y = lift((t, v) -> matched_trace_vectors(t, v, offsets["presynaptic_si3"])[2], trace_obs[("presynaptic", :time)], trace_obs[("presynaptic", :V3)])
    post_si3_x = lift((t, v) -> matched_trace_vectors(t, v, offsets["postsynaptic_si3"])[1], trace_obs[("postsynaptic", :time)], trace_obs[("postsynaptic", :V3)])
    post_si3_y = lift((t, v) -> matched_trace_vectors(t, v, offsets["postsynaptic_si3"])[2], trace_obs[("postsynaptic", :time)], trace_obs[("postsynaptic", :V3)])
    lines!(axE, pre_drive_x, pre_drive_y; color = RGBf(0.35, 0.35, 0.35), linewidth = 1.4, visible = scan_protocol_visible)
    lines!(axE, pre_si2_x, pre_si2_y; color = DASH_COLORS["presynaptic"], linewidth = 1.4)
    lines!(axE, post_si2_x, post_si2_y; color = DASH_COLORS["postsynaptic"], linewidth = 1.4)
    lines!(axE, shared_sm_x, shared_sm_y; color = RGBf(0.45, 0.24, 0.62), linewidth = 1.1)
    lines!(axE, pre_si3_x, pre_si3_y; color = DASH_COLORS["presynaptic"], linewidth = 1.4, visible = cb_show_si3.checked)
    lines!(axE, post_si3_x, post_si3_y; color = DASH_COLORS["postsynaptic"], linewidth = 1.4, visible = cb_show_si3.checked)
    for mode_name in MODE_LABELS
        si3_burst_markers_visible = lift((show_si3, show_markers) -> show_si3 && show_markers, cb_show_si3.checked, cb_show_burst_markers.checked)
        for (signal, offset_key, visible) in ((:V1, mode_name, cb_show_burst_markers.checked), (:V3, "$(mode_name)_si3", si3_burst_markers_visible))
            time_obs = trace_obs[(mode_name, :time)]
            voltage_obs = trace_obs[(mode_name, signal)]
            tick_x = lift((t, v, d) -> burst_tick_vectors(t, v, offsets[offset_key]; detector_kwargs(d)...)[1], time_obs, voltage_obs, detector_obs)
            tick_y = lift((t, v, d) -> burst_tick_vectors(t, v, offsets[offset_key]; detector_kwargs(d)...)[2], time_obs, voltage_obs, detector_obs)
            marker_x = lift((t, v, d) -> burst_marker_vectors(t, v, offsets[offset_key]; detector_kwargs(d)...)[1], time_obs, voltage_obs, detector_obs)
            marker_y = lift((t, v, d) -> burst_marker_vectors(t, v, offsets[offset_key]; detector_kwargs(d)...)[2], time_obs, voltage_obs, detector_obs)
            lines!(axE, tick_x, tick_y; color = (DASH_COLORS[mode_name], 0.28), linewidth = 1.0, visible = visible)
            scatter!(axE, marker_x, marker_y; color = DASH_COLORS[mode_name], marker = DASH_MARKERS[mode_name], markersize = 10, strokecolor = :black, strokewidth = 0.8, visible = visible)
        end
    end
    x_min, x_max, y_min, y_max = trace_panel_limits(bio_t1_obs[], bio_v1_obs[], bio_t2_obs[], bio_v2_obs[], traces_obs[]; show_si3 = cb_show_si3.checked[])
    label_x = x_min + 0.045 * (x_max - x_min)
    text!(axE, label_x, offsets["si1"] + 38; text = "Si1 biology", align = (:left, :center), fontsize = 15, color = RGBf(0.35, 0.35, 0.35), visible = empirical_protocol_visible)
    text!(axE, label_x, offsets["si2"] + 38; text = "Si2 biology", align = (:left, :center), fontsize = 15, color = RGBf(0.1, 0.1, 0.1), visible = empirical_protocol_visible)
    text!(axE, label_x, offsets["si1"] + 38; text = "Si1 scan protocol", align = (:left, :center), fontsize = 15, color = RGBf(0.35, 0.35, 0.35), visible = scan_protocol_visible)
    text!(axE, label_x, offsets["presynaptic"] + 38; text = "Si2 presynaptic", align = (:left, :center), fontsize = 15, color = DASH_COLORS["presynaptic"])
    text!(axE, label_x, offsets["presynaptic_si3"] + 38; text = "Si3 presynaptic", align = (:left, :center), fontsize = 15, color = DASH_COLORS["presynaptic"], visible = cb_show_si3.checked)
    text!(axE, label_x, offsets["postsynaptic"] + 38; text = "Si2 postsynaptic", align = (:left, :center), fontsize = 15, color = DASH_COLORS["postsynaptic"])
    text!(axE, label_x, shared_neuromod_offset(offsets) + 22; text = "sm neuromodulation", align = (:left, :center), fontsize = 12, color = RGBf(0.45, 0.24, 0.62))
    text!(axE, label_x, offsets["postsynaptic_si3"] + 38; text = "Si3 postsynaptic", align = (:left, :center), fontsize = 15, color = DASH_COLORS["postsynaptic"], visible = cb_show_si3.checked)
    trace_panel_title(source, selected, calcium_ics) = begin
        gain = round(selected[2], digits = 3)
        protocol = source == :scan_protocol ? "scan burst protocol" : "empirical Si1"
        ca = "Ca IC=$(round(calcium_ics.Ca0, digits = 2))/$(round(calcium_ics.Ca1, digits = 2))/$(round(calcium_ics.Ca2, digits = 2))/$(round(calcium_ics.Ca3, digits = 2))/$(round(calcium_ics.Ca4, digits = 2))"
        "Representative Traces ($(selected[1]), gain=$(gain), $(protocol), $(ca))"
    end
    axE.ylabel = "V (mV)"
    axE.xticksvisible = false
    axE.yticksvisible = false
    axE.xticklabelsvisible = false
    axE.yticklabelsvisible = false
    axE.xgridvisible = false
    axE.ygridvisible = false
    hidespines!(axE, :t, :r)
    add_l_scale_bar!(axE, x_min, x_max, y_min, y_max)
    reset_trace_limits!(axE, bio_t1_obs[], bio_v1_obs[], bio_t2_obs[], bio_v2_obs[], traces_obs[]; show_si3 = cb_show_si3.checked[])
    add_panel_label!(axE, "E")

    param_value(key::Symbol) = Float64(sg_params.by_key[key].value[])
    function current_param_values(mode::ControlMode)
        dashboard_params(
            mode;
            x_shift_si2 = param_value(:x_shift_si2),
            x_shift_si3 = param_value(:x_shift_si3),
            presyn_g0 = param_value(:presyn_g0),
            postsyn_g0 = param_value(:postsyn_g0),
            alpha1 = param_value(:alpha1),
            beta1 = param_value(:beta1),
            s0_floor = param_value(:s0_floor),
            alpham = param_value(:alpham),
            betam = param_value(:betam),
            sm_floor = param_value(:sm_floor),
            ca_shift_si2 = param_value(:ca_shift_si2),
            ca_shift_si3 = param_value(:ca_shift_si3),
            presynaptic_base_g = param_value(:presynaptic_base_g),
            direct_post_base_g = param_value(:direct_post_base_g),
            presyn_alphax = param_value(:presyn_alphax),
            postsyn_alphax = param_value(:postsyn_alphax),
            presyn_betax = param_value(:presyn_betax),
            postsyn_betax = param_value(:postsyn_betax),
            si3_exc_floor = param_value(:si3_exc_floor),
            slow_inhib_g = param_value(:slow_inhib_g),
            slow_inhib_alpha = param_value(:slow_inhib_alpha),
            slow_inhib_beta = param_value(:slow_inhib_beta),
            si2_mutual_inhib_g = param_value(:si2_mutual_inhib_g),
            si2_mutual_inhib_alpha = param_value(:si2_mutual_inhib_alpha),
            si2_mutual_inhib_beta = param_value(:si2_mutual_inhib_beta),
            si_inhib_floor = param_value(:si_inhib_floor),
            si3_mutual_inhib_g = param_value(:si3_mutual_inhib_g),
            si3_mutual_inhib_alpha = param_value(:si3_mutual_inhib_alpha),
            si3_mutual_inhib_beta = param_value(:si3_mutual_inhib_beta),
            t1_ms = param_value(:t1_ms),
        )
    end
    current_config() = dashboard_config(Float64(sg_scan.sliders[4].value[]), Float64(sg_scan.sliders[5].value[]); detector_kwargs(current_burst_detector())...)
    current_trace_source() = cb_trace_scan_protocol.checked[] ? :scan_protocol : :empirical
    current_trace_time_scale() = Float64(sg_timescale.sliders[1].value[])
    current_trace_calcium_ics() = (
        Ca0 = Float64(sg_trace.sliders[3].value[]),
        Ca1 = Float64(sg_trace.sliders[4].value[]),
        Ca2 = Float64(sg_trace.sliders[5].value[]),
        Ca3 = Float64(sg_trace.sliders[6].value[]),
        Ca4 = Float64(sg_trace.sliders[7].value[]),
    )
    current_burst_detector() = (
        spike_threshold_mv = Float64(sg_trace.sliders[8].value[]),
        spike_refractory_s = Float64(sg_trace.sliders[9].value[]),
        burst_factor = Float64(sg_trace.sliders[10].value[]),
    )
    current_param_default_values() = [param_value(key) for key in DASHBOARD_PARAM_KEYS]

    function selected_terminal_state()
        selected = selected_obs[]
        trace = get(traces_obs[], selected[1], DataFrame())
        isempty(trace) && return nothing
        all(col -> col in propertynames(trace), state_columns()) || return nothing
        last_idx = nrow(trace)
        terminal = [Float64(trace[last_idx, col]) for col in state_columns()]
        any(!isfinite, terminal) && return nothing
        return terminal
    end

    function sync_calcium_sliders_to_state!(u::Vector{Float64})
        length(u) >= 15 || return nothing
        for (slider, idx) in zip(sg_trace.sliders[3:7], 11:15)
            set_close_to!(slider, clamp(Float64(u[idx]), first(slider.range[]), last(slider.range[])))
        end
        return nothing
    end

    function reset_current_trace_limits!()
        source = current_trace_source()
        reset_trace_limits!(axE, bio_t1_obs[], bio_v1_obs[], bio_t2_obs[], bio_v2_obs[], traces_obs[]; source = source == :scan_protocol ? :scan_protocol : :representative, show_si3 = cb_show_si3.checked[])
        return nothing
    end

    function refresh_time_scaled_biology!()
        time_scale = current_trace_time_scale()
        trace_time_scale_obs[] = time_scale
        t1, v1, t2, v2 = load_biology_pair(; time_scale = time_scale)
        bio_t1_obs[] = t1
        bio_v1_obs[] = v1
        bio_t2_obs[] = t2
        bio_v2_obs[] = v2
        lx, ly = load_biological_fig2c(; time_scale = time_scale)
        lit_x_obs[] = lx
        lit_y_obs[] = ly
        bio_points_obs[] = load_biological_trace_points(; time_scale = time_scale, detector_kwargs(current_burst_detector())...)
        reset_current_trace_limits!()
        return nothing
    end

    function refresh_detector_outputs!()
        detector = current_burst_detector()
        detector_obs[] = detector
        bio_points_obs[] = load_biological_trace_points(; time_scale = current_trace_time_scale(), detector_kwargs(detector)...)
        model_points_obs[] = representative_model_points(traces_obs[]; detector_kwargs(detector)...)
        refresh_model_point_observables!(model_points_obs, point_obs)
        if current_trace_source() == :scan_protocol
            request_trace_update!(:all)
        end
        status_obs[] = "Burst detector: spike threshold=$(round(detector.spike_threshold_mv, digits = 1)), refractory=$(round(detector.spike_refractory_s, digits = 3))s, factor=$(round(detector.burst_factor, digits = 2))"
        return nothing
    end
    axE.title[] = trace_panel_title(current_trace_source(), selected_obs[], current_trace_calcium_ics())

    function run_trace_for_mode!(mode::ControlMode, gain::Float64)
        mode_name = mode_label(mode)
        source = current_trace_source()
        source_label = source == :scan_protocol ? "scan burst protocol" : "empirical Si1"
        status_obs[] = "Running representative trace ($(source_label)): $(mode_name) gain=$(round(gain, digits = 4))"
        calcium_ics = current_trace_calcium_ics()
        trace = if source == :scan_protocol
            simulate_scan_state_trace(
                current_param_values(mode),
                gain,
                mode,
                current_config();
                display_time_scale = current_trace_time_scale(),
                calcium_ics = calcium_ics,
                initial_u0 = full_initial_state[],
            )
        else
            simulate_driven_state_trace(
                current_param_values(mode),
                gain,
                mode,
                bio_t1_obs[],
                bio_v1_obs[];
                saveat_ms = Float64(sg_scan.sliders[5].value[]),
                display_time_scale = 1.0,
                calcium_ics = calcium_ics,
                initial_u0 = full_initial_state[],
            )
        end
        traces = copy(traces_obs[])
        traces[mode_name] = trace
        traces_obs[] = traces
        set_model_trace_observables!(trace_obs, traces)
        model_points_obs[] = representative_model_points(traces; detector_kwargs(current_burst_detector())...)
        refresh_model_point_observables!(model_points_obs, point_obs)
        reset_current_trace_limits!()
        status_obs[] = "Trace complete: $(mode_name) gain=$(round(gain, digits = 4))"
        return trace
    end

    function run_selected_trace!()
        selected = selected_obs[]
        run_trace_for_mode!(parse_mode_label(selected[1]), current_gain_for_mode(selected[1]))
        return traces_obs[]
    end

    function current_gain_for_mode(mode_name::String)
        mode_name == "presynaptic" && return Float64(sg_trace.sliders[1].value[])
        mode_name == "postsynaptic" && return Float64(sg_trace.sliders[2].value[])
        error("Unknown representative trace mode $(mode_name)")
    end

    compact_param(x; digits = 4) = string(round(Float64(x), digits = digits))

    function active_mode_summary(mode::ControlMode, gain::Float64)
        params = current_param_values(mode)
        sim_params = mode_params(params, mode)
        active_alpha = baseline_alphax(sim_params, mode)
        active_g41 = baseline_g41(sim_params, mode)
        active_g32 = baseline_g32(sim_params, mode)
        return "$(mode_label(mode)): gain=$(compact_param(gain, digits = 3)), g0=$(compact_param(sim_params.g0)), " *
               "alpha_x=$(compact_param(active_alpha, digits = 3)), g41=$(compact_param(active_g41)), " *
               "g32=$(compact_param(active_g32)), beta_x=$(compact_param(sim_params.betax)), " *
               "Ca1/3=$(compact_param(sim_params.Ca_shift1, digits = 1))/$(compact_param(sim_params.Ca_shift3, digits = 1))"
    end

    function refresh_active_parameter_readout!()
        active_param_text[] = "Active baselines\n" *
            active_mode_summary(presynaptic, current_gain_for_mode("presynaptic")) * "\n" *
            active_mode_summary(postsynaptic, current_gain_for_mode("postsynaptic"))
        return nothing
    end

    refresh_active_parameter_readout!()

    function run_all_traces!()
        run_trace_for_mode!(presynaptic, current_gain_for_mode("presynaptic"))
        run_trace_for_mode!(postsynaptic, current_gain_for_mode("postsynaptic"))
        return traces_obs[]
    end

    trace_dirty = Ref(false)
    trace_running = Ref(false)
    trace_scope = Ref(:selected)
    function request_trace_update!(scope::Symbol = :selected)
        if scope == :all || trace_scope[] != :all
            trace_scope[] = scope
        end
        trace_dirty[] = true
        trace_running[] && return
        @async begin
            trace_running[] = true
            try
                while trace_dirty[]
                    trace_dirty[] = false
                    scope = trace_scope[]
                    trace_scope[] = :selected
                    try
                        if scope == :all
                            run_all_traces!()
                        elseif scope == :presynaptic
                            run_trace_for_mode!(presynaptic, current_gain_for_mode("presynaptic"))
                        elseif scope == :postsynaptic
                            run_trace_for_mode!(postsynaptic, current_gain_for_mode("postsynaptic"))
                        else
                            run_selected_trace!()
                        end
                    catch err
                        status_obs[] = "Trace failed: $(sprint(showerror, err))"
                        rethrow()
                    end
                    yield()
                end
            finally
                trace_running[] = false
            end
        end
        return nothing
    end

    function reset_params!()
        for (key, value) in zip(DASHBOARD_PARAM_KEYS, param_defaults)
            set_close_to!(sg_params.by_key[key], value)
        end
        refresh_active_parameter_readout!()
        request_trace_update!(:all)
        status_obs[] = "Parameters reset"
        return nothing
    end

    function set_initial_conditions_to_terminal!()
        terminal = selected_terminal_state()
        if isnothing(terminal)
            status_obs[] = "No complete terminal state available for selected trace"
            return nothing
        end
        full_initial_state[] = terminal
        sync_calcium_sliders_to_state!(terminal)
        axE.title[] = trace_panel_title(current_trace_source(), selected_obs[], current_trace_calcium_ics())
        request_trace_update!(:all)
        status_obs[] = "Set full initial state from selected terminal state"
        return nothing
    end

    function match_post_baseline_to_pre!()
        set_close_to!(sg_params.by_key[:postsyn_g0], param_value(:presyn_g0))
        set_close_to!(sg_params.by_key[:postsyn_alphax], param_value(:presyn_alphax))
        set_close_to!(sg_params.by_key[:postsyn_betax], param_value(:presyn_betax))
        set_close_to!(sg_params.by_key[:direct_post_base_g], param_value(:presynaptic_base_g))
        refresh_active_parameter_readout!()
        request_trace_update!(:all)
        status_obs[] = "Matched postsynaptic baseline g0/alpha_x/g41/g32 to presynaptic"
        return nothing
    end

    function run_selected_scan!()
        modes = ControlMode[]
        cb_presyn_scan.checked[] && push!(modes, presynaptic)
        cb_postsyn_scan.checked[] && push!(modes, postsynaptic)
        isempty(modes) && (status_obs[] = "No scan mode selected"; return)
        n = Int(round(sg_scan.sliders[3].value[]))
        gains = ordered_gain_grid(Float64(sg_scan.sliders[1].value[]), Float64(sg_scan.sliders[2].value[]), n, String(order_menu.selection[]))
        config = current_config()
        keep_modes = Set(mode_label.(modes))
        raw_acc = copy(raw_obs[][[!(m in keep_modes) for m in raw_obs[].mode], :])
        summary_acc = copy(summary_obs[][[!(m in keep_modes) for m in summary_obs[].mode], :])
        total = length(modes) * length(gains)
        completed = 0
        for mode in modes
            params = current_param_values(mode)
            for gain in gains
                completed += 1
                status_obs[] = "Scanning $(mode_label(mode)) gain=$(round(gain, digits = 4)) ($(completed)/$(total))"
                raw, summary = run_scan_point(params, mode, gain, config)
                append!(raw_acc, raw; cols = :union)
                append!(summary_acc, summary; cols = :union)
                raw_obs[] = copy(raw_acc)
                summary_obs[] = copy(summary_acc)
                refresh_scan_observables!(summary_obs, raw_obs, plot_obs)
                yield()
            end
        end
        status_obs[] = "Scan complete ($(total) simulations)"
    end

    on(run_scan_button.clicks) do _
        @async run_selected_scan!()
    end
    on(reset_params_button.clicks) do _
        reset_params!()
    end
    on(match_baselines_button.clicks) do _
        match_post_baseline_to_pre!()
    end
    on(terminal_ic_button.clicks) do _
        set_initial_conditions_to_terminal!()
    end
    on(save_params_button.clicks) do _
        path = save_param_defaults(current_param_default_values())
        param_defaults = load_param_defaults()
        status_obs[] = "Saved parameter defaults to $(path)"
    end
    on(save_button.clicks) do _
        trace_source = current_trace_source() == :scan_protocol ? :scan_protocol : :representative
        clean = build_clean_figure(summary_obs[], raw_obs[], traces_obs[]; trace_source = trace_source, show_si3 = cb_show_si3.checked[], show_burst_markers = cb_show_burst_markers.checked[], detector = current_burst_detector(), time_scale = current_trace_time_scale())
        mkpath(DASHBOARD_OUTPUT_DIR)
        CairoMakie.save(DASHBOARD_CLEAN_PNG, clean)
        status_obs[] = "Saved $(DASHBOARD_CLEAN_PNG)"
    end
    on(write_csv_button.clicks) do _
        mkpath(DASHBOARD_OUTPUT_DIR)
        CSV.write(DASHBOARD_POINTS_CSV, raw_obs[])
        CSV.write(DASHBOARD_SUMMARY_CSV, summary_obs[])
        status_obs[] = "Wrote dashboard CSVs"
    end
    on(sg_trace.sliders[1].value) do gain
        selected = selected_obs[]
        selected[1] == "presynaptic" && (selected_obs[] = ("presynaptic", Float64(gain)))
        refresh_active_parameter_readout!()
        request_trace_update!(:presynaptic)
    end
    on(sg_trace.sliders[2].value) do gain
        selected = selected_obs[]
        selected[1] == "postsynaptic" && (selected_obs[] = ("postsynaptic", Float64(gain)))
        refresh_active_parameter_readout!()
        request_trace_update!(:postsynaptic)
    end
    for slider in sg_trace.sliders[3:7]
        on(slider.value) do _
            axE.title[] = trace_panel_title(current_trace_source(), selected_obs[], current_trace_calcium_ics())
            request_trace_update!(:all)
        end
    end
    for slider in sg_trace.sliders[8:10]
        on(slider.value) do _
            refresh_detector_outputs!()
        end
    end
    on(sg_timescale.sliders[1].value) do value
        refresh_time_scaled_biology!()
        axE.title[] = trace_panel_title(current_trace_source(), selected_obs[], current_trace_calcium_ics())
        request_trace_update!(:all)
        status_obs[] = "Trace time scale=$(round(Float64(value), digits = 2))x"
    end
    for slider in sg_params.sliders
        on(slider.value) do _
            refresh_active_parameter_readout!()
            request_trace_update!(:all)
        end
    end
    on(sg_scan.sliders[5].value) do _
        request_trace_update!(:all)
    end
    on(cb_trace_scan_protocol.checked) do _
        source = current_trace_source()
        axE.title[] = trace_panel_title(source, selected_obs[], current_trace_calcium_ics())
        reset_current_trace_limits!()
        request_trace_update!(:all)
    end
    on(cb_show_si3.checked) do show_si3
        reset_current_trace_limits!()
    end
    on(selected_obs) do selected
        axE.title[] = trace_panel_title(current_trace_source(), selected, current_trace_calcium_ics())
    end
    on(events(fig).mousebutton) do event
        event.button == Mouse.left || return
        event.action == Mouse.press || return
        plt, idx = pick(fig)
        idx < 1 && return
        haskey(pick_targets, plt) || return
        mode_name, _ = pick_targets[plt]
        sub = sort(summary_obs[][summary_obs[].mode .== mode_name, :], :control_gain)
        idx > nrow(sub) && return
        gain = Float64(sub[idx, :control_gain])
        selected_obs[] = (mode_name, gain)
        slider = mode_name == "presynaptic" ? sg_trace.sliders[1] : sg_trace.sliders[2]
        set_close_to!(slider, gain)
        request_trace_update!(Symbol(mode_name))
    end

    if cb_show_si3.checked[] && !traces_have_finite_si3(traces_obs[])
        status_obs[] = "Startup cache lacks Si3 traces; recalculating representative traces"
        request_trace_update!(:all)
    end

    colsize!(fig.layout, 1, Fixed(360))
    colsize!(fig.layout, 2, Fixed(410))
    rowsize!(plots, 1, Fixed(240))
    rowsize!(plots, 2, Fixed(330))
    rowsize!(plots, 3, Fixed(560))
    rowgap!(fig.layout, 14)
    colgap!(fig.layout, 14)

    handle = (
        fig = fig,
        raw = raw_obs,
        summary = summary_obs,
        selected = selected_obs,
        traces = traces_obs,
        model_points = model_points_obs,
        bio_points = bio_points_obs,
        status = status_obs,
        plot_obs = plot_obs,
        trace_obs = trace_obs,
        scan_controls = sg_scan,
        trace_controls = sg_trace,
        time_scale_control = sg_timescale,
        param_controls = sg_params,
        order_menu = order_menu,
        checkboxes = (
            scan_presynaptic = cb_presyn_scan,
            scan_postsynaptic = cb_postsyn_scan,
            trace_scan_protocol = cb_trace_scan_protocol,
            show_si3 = cb_show_si3,
            show_burst_markers = cb_show_burst_markers,
        ),
        buttons = (
            run_scan = run_scan_button,
            set_terminal_ics = terminal_ic_button,
            match_baselines = match_baselines_button,
            reset_params = reset_params_button,
            save_params = save_params_button,
            save = save_button,
            write_csv = write_csv_button,
        ),
        axes = (circuit = axA, compare = axB, during = axC, after = axD, traces = axE),
        actions = (
            run_scan! = run_selected_scan!,
            run_trace! = run_selected_trace!,
            run_all_traces! = run_all_traces!,
            run_trace_for_mode! = run_trace_for_mode!,
            request_trace_update! = request_trace_update!,
            reset_params! = reset_params!,
            current_config = current_config,
            current_params = current_param_values,
            current_trace_source = current_trace_source,
            current_trace_time_scale = current_trace_time_scale,
            current_trace_calcium_ics = current_trace_calcium_ics,
            selected_terminal_state = selected_terminal_state,
            set_initial_conditions_to_terminal! = set_initial_conditions_to_terminal!,
            current_burst_detector = current_burst_detector,
            current_param_default_values = current_param_default_values,
            refresh_detector_outputs! = refresh_detector_outputs!,
            refresh_time_scaled_biology! = refresh_time_scaled_biology!,
            reset_trace_limits! = reset_current_trace_limits!,
        ),
    )
    DASHBOARD_HANDLE[] = handle
    return handle
end
