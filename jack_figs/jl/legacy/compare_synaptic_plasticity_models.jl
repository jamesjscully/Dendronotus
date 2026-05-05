using CSV
using DataFrames
using OrdinaryDiffEq
using Plots
using Plots.PlotMeasures
using Statistics

gr()

const SCRIPT_DIR = @__DIR__
const OUTPUT_DIR = joinpath(dirname(SCRIPT_DIR), "outputs")
const RAW_CSV = joinpath(OUTPUT_DIR, "plasticity_comparison_points.csv")
const SUMMARY_CSV = joinpath(OUTPUT_DIR, "plasticity_comparison_summary.csv")
const FIGURE_PNG = joinpath(OUTPUT_DIR, "plasticity_comparison_scatter.png")

@enum ControlMode direct_gmax presynaptic postsynaptic

struct SweepConfig
    max_time_s::Float64
    saveat_ms::Float64
    spike_threshold_mv::Float64
    spike_refractory_s::Float64
    burst_factor::Float64
    onset_timeout_s::Float64
    transient_bursts::Int
    measured_cycles::Int
end

struct ModeSweepSpec
    mode::ControlMode
    label::String
    control_values::Vector{Float64}
    control_axis_label::String
end

struct ModelParams
    direct_post_base_g::Float64
    presynaptic_base_alpha::Float64
    presynaptic_base_g::Float64
    gelec::Float64
    alpha1::Float64
    beta1::Float64
    scale1::Float64
    s0_floor::Float64
    g0::Float64
    alpham::Float64
    betam::Float64
    scale_sm::Float64
    sm_floor::Float64
    alphax::Float64
    betax::Float64
    g41::Float64
    g32::Float64
    g14::Float64
    g23::Float64
    alphai::Float64
    betai::Float64
    scalei::Float64
    si_inhib_floor::Float64
    alpha3::Float64
    beta3::Float64
    scale3::Float64
    g34::Float64
    g43::Float64
    si3_exc_floor::Float64
    alpha2::Float64
    beta2::Float64
    scale2::Float64
    g21::Float64
    g12::Float64
    Ca_shift0::Float64
    Ca_shift1::Float64
    Ca_shift2::Float64
    Ca_shift3::Float64
    Ca_shift4::Float64
    x_shift::Float64
    x_shift1::Float64
    x_shift2::Float64
    x_shift3::Float64
    x_shift4::Float64
    Iapp::Float64
    t1_ms::Float64
    t2_ms::Float64
end

function default_config()
    SweepConfig(
        9000.0,
        2.0,
        -20.0,
        0.02,
        2.0,
        200.0,
        3,
        10,
    )
end

mode_label(mode::ControlMode) = string(mode)
function parse_mode_label(label::AbstractString)
    label == "direct_gmax" && return direct_gmax
    label == "presynaptic" && return presynaptic
    label == "postsynaptic" && return postsynaptic
    error("Unknown mode label `$(label)`")
end

function default_sweep_specs()
    [
        ModeSweepSpec(
            direct_gmax,
            "Direct g_max",
            collect(range(0.0, 6.0; length = 50)),
            "Direct g_max Gain",
        ),
        ModeSweepSpec(
            presynaptic,
            "Presynaptic alpha",
            collect(range(0.0, 6.0; length = 60)),
            "Presynaptic alpha Gain",
        ),
        ModeSweepSpec(
            postsynaptic,
            "Postsynaptic g_syn",
            collect(range(0.0, 6.0; length = 50)),
            "Postsynaptic g_syn Gain",
        ),
    ]
end

function default_params()
    alpha1 = 0.01
    beta1 = 0.002
    s0_floor = 0.01
    alpham = 0.005
    betam = 0.0001
    sm_floor = 0.01
    alphax = 0.05
    betax = 0.001
    alphai = 0.012
    betai = 0.001
    si_inhib_floor = 0.01
    alpha3 = 0.01
    beta3 = 0.004
    si3_exc_floor = 0.01
    alpha2 = 0.01
    beta2 = 0.005

    ModelParams(
        0.000625,
        0.011,
        0.00125,
        0.002,
        alpha1,
        beta1,
        alpha1 / (alpha1 + beta1),
        s0_floor,
        0.002,
        alpham,
        betam,
        (alpham - betam) / alpham,
        sm_floor,
        alphax,
        betax,
        0.000625,
        0.000625,
        0.01,
        0.01,
        alphai,
        betai,
        (alphai - betai) / alphai,
        si_inhib_floor,
        alpha3,
        beta3,
        alpha3 / (alpha3 + beta3),
        0.02,
        0.02,
        si3_exc_floor,
        alpha2,
        beta2,
        alpha2 / (alpha2 + beta2),
        0.02,
        0.02,
        -10.0,
        -25.0,
        -25.0,
        -25.0,
        -25.0,
        -4.0,
        -4.0,
        -4.0,
        -4.0,
        -4.0,
        0.1,
        1000.0,
        35000.0,
    )
end

function updated_params(params::ModelParams; kwargs...)
    names = fieldnames(ModelParams)
    lookup = Dict(kwargs)
    values = Any[]
    for name in names
        push!(values, get(lookup, name, getfield(params, name)))
    end
    return ModelParams(values...)
end

function mode_params(params::ModelParams, mode::ControlMode)
    if mode == presynaptic
        return updated_params(
            params;
            presynaptic_base_alpha = params.presynaptic_base_alpha,
            presynaptic_base_g = params.presynaptic_base_g,
            g0 = params.g0,
            alpha3 = params.alpha3,
            beta3 = params.beta3,
            scale3 = params.scale3,
            g34 = params.g34,
            g43 = params.g43,
            alpha2 = params.alpha2,
            beta2 = params.beta2,
            scale2 = params.scale2,
            g21 = params.g21,
            g12 = params.g12,
        )
    end

    return params
end

heaviside(x) = x >= 0 ? 1.0 : 0.0
syn_activation(v) = 1.0 / (1.0 + exp(-20.0 * (v + 20.0)))

baseline_alphax(params::ModelParams, mode::ControlMode) =
    mode == presynaptic ? params.presynaptic_base_alpha : params.alphax

baseline_g41(params::ModelParams, mode::ControlMode) =
    mode == presynaptic ? params.presynaptic_base_g : params.direct_post_base_g

baseline_g32(params::ModelParams, mode::ControlMode) =
    mode == presynaptic ? params.presynaptic_base_g : params.direct_post_base_g

function initial_state()
    u0 = zeros(37)
    u0[1:5] .= -44.0
    u0[6:10] .= [0.9, 0.6, 0.6, 0.5, 0.6]
    u0[11:15] .= [0.3, 1.0, 1.0, 1.1, 1.0]
    u0[16:25] .= 0.0
    u0[26:27] .= 0.0
    u0[28:32] .= 0.0
    u0[33:36] .= [0.0, 0.0, 0.1, 0.1]
    u0[37] = 0.5
    return u0
end

mutable struct SpikeTracker
    spike_times_ms::Vector{Float64}
    last_spike_ms::Float64
    spike_refractory_ms::Float64
    threshold_mv::Float64
    required_complete_bursts::Int
end

mutable struct ProtocolTracker
    phase::Symbol
    current_on::Bool
    onset_failed::Bool
    stimulus_off_s::Float64
    withdrawal_end_s::Float64
    measured_start_s::Float64
    measured_stop_s::Float64
end

struct BurstDetection
    burst_factor::Float64
    burst_isi_indices::Vector{Int}
    begin_indices::Vector{Int}
    end_indices::Vector{Int}
    start_times_s::Vector{Float64}
    end_times_s::Vector{Float64}
    spike_counts::Vector{Int}
    is_complete::Vector{Bool}
end

function make_spike_callback(tracker::SpikeTracker)
    condition(u, t, integrator) = u[2] - tracker.threshold_mv

    function affect!(integrator)
        t_ms = integrator.t
        if isempty(tracker.spike_times_ms) || (t_ms - tracker.last_spike_ms) >= tracker.spike_refractory_ms
            push!(tracker.spike_times_ms, t_ms)
            tracker.last_spike_ms = t_ms
            if tracker.required_complete_bursts > 0
                bursts = detect_bursts_from_spikes(tracker.spike_times_ms ./ 1000.0, integrator.p.config)
                if count(identity, bursts.is_complete) >= tracker.required_complete_bursts
                    terminate!(integrator)
                end
            end
        end
    end

    return ContinuousCallback(condition, affect!; affect_neg! = nothing, save_positions = (false, false))
end

function make_protocol_callback(protocol::ProtocolTracker, terminate_on_withdrawal_end::Bool = true)
    condition(u, t, integrator) = begin
        if protocol.phase == :seeking_onset
            return t - (integrator.p.params.t1_ms + integrator.p.config.onset_timeout_s * 1000.0)
        elseif terminate_on_withdrawal_end && protocol.phase == :withdrawal
            return t - protocol.withdrawal_end_s * 1000.0
        else
            return 1.0
        end
    end

    function affect!(integrator)
        if protocol.phase == :seeking_onset
            protocol.onset_failed = true
            protocol.current_on = false
            protocol.phase = :failed
            terminate!(integrator)
        elseif protocol.phase == :withdrawal
            protocol.phase = :done
            terminate!(integrator)
        end
    end

    return ContinuousCallback(condition, affect!; affect_neg! = nothing, save_positions = (false, false))
end

function update_protocol_from_spikes!(protocol::ProtocolTracker, spike_times_ms::Vector{Float64}, config::SweepConfig, integrator)
    bursts = detect_bursts_from_spikes(spike_times_ms ./ 1000.0, config)
    complete_count = count(identity, bursts.is_complete)

    if protocol.phase == :seeking_onset && complete_count >= 1
        protocol.phase = :stimulated
    end

    if protocol.phase != :withdrawal && protocol.phase != :done && complete_count >= config.transient_bursts + config.measured_cycles
        complete_starts = complete_burst_start_times(bursts)
        complete_ends = complete_burst_end_times(bursts)
        measured_start_idx = config.transient_bursts + 1
        measured_stop_idx = config.transient_bursts + config.measured_cycles
        protocol.measured_start_s = complete_starts[measured_start_idx]
        protocol.measured_stop_s = complete_ends[measured_stop_idx]
        protocol.stimulus_off_s = integrator.t / 1000.0
        withdrawal_duration_s = max(protocol.measured_stop_s - protocol.measured_start_s, 0.0)
        protocol.withdrawal_end_s = protocol.stimulus_off_s + withdrawal_duration_s
        protocol.current_on = false
        protocol.phase = :withdrawal
        add_tstop!(integrator, protocol.withdrawal_end_s * 1000.0)
    end

    return nothing
end

function empty_burst_detection()
    BurstDetection(NaN, Int[], Int[], Int[], Float64[], Float64[], Int[], Bool[])
end

function detect_bursts_from_spikes(spike_times_s::AbstractVector{<:Real}, config::SweepConfig)
    times = sort(Float64.(collect(spike_times_s)))
    length(times) < 4 && return empty_burst_detection()

    isi_values = diff(times)
    length(isi_values) < 3 && return empty_burst_detection()

    burst_isi_indices = Int[]
    last_cut_idx = -1
    for i in 2:(length(isi_values) - 1)
        start_idx = last_cut_idx + 2
        end_idx = i - 1
        end_idx < start_idx && continue

        isi_window = isi_values[start_idx:end_idx]
        isempty(isi_window) && continue
        d_median = median(isi_window)

        if isi_values[i] > (config.burst_factor * d_median) &&
           isi_values[i + 1] < (isi_values[i] / config.burst_factor)
            push!(burst_isi_indices, i + 1)
            last_cut_idx = i - 1
        end
    end

    isempty(burst_isi_indices) && return empty_burst_detection()

    start_times_s = Float64[]
    end_times_s = Float64[]
    spike_counts = Int[]
    is_complete = Bool[]
    kept_begin_indices = Int[]
    kept_end_indices = Int[]
    boundaries = vcat(burst_isi_indices, length(times))

    first_end_idx = boundaries[1]
    if first_end_idx > 1
        push!(kept_begin_indices, 1)
        push!(kept_end_indices, first_end_idx)
        push!(start_times_s, times[1])
        push!(end_times_s, times[first_end_idx])
        push!(spike_counts, first_end_idx)
        push!(is_complete, true)
    end

    for boundary_idx in 1:(length(boundaries) - 1)
        prev_idx = boundaries[boundary_idx]
        next_idx = boundaries[boundary_idx + 1]
        next_idx - prev_idx <= 1 && continue
        begin_idx = prev_idx + 1
        end_idx = next_idx
        push!(kept_begin_indices, begin_idx)
        push!(kept_end_indices, end_idx)
        push!(start_times_s, times[begin_idx])
        push!(end_times_s, times[end_idx])
        push!(spike_counts, end_idx - begin_idx + 1)
        push!(is_complete, boundary_idx < (length(boundaries) - 1))
    end

    isempty(start_times_s) && return empty_burst_detection()
    count(identity, is_complete) < 2 && return empty_burst_detection()
    return BurstDetection(config.burst_factor, burst_isi_indices, kept_begin_indices, kept_end_indices, start_times_s, end_times_s, spike_counts, is_complete)
end

complete_burst_start_times(bursts::BurstDetection) = bursts.start_times_s[bursts.is_complete]
complete_burst_end_times(bursts::BurstDetection) = bursts.end_times_s[bursts.is_complete]
incomplete_burst_start_times(bursts::BurstDetection) = bursts.start_times_s[.!bursts.is_complete]

function network_ode!(du, u, p, t_ms)
    params = p.params
    control_gain = p.control_gain
    mode = p.mode

    V0, V1, V2, V3, V4 = u[1:5]
    x0, x1, x2, x3, x4 = u[6:10]
    Ca0, Ca1, Ca2, Ca3, Ca4 = u[11:15]
    h0, h1, h2, h3, h4 = u[16:20]
    n0, n1, n2, n3, n4 = u[21:25]
    y1, y2 = u[26:27]
    s0, s1, s2, s3, s4 = u[28:32]
    s12, s21, s34, s43 = u[33:36]
    sm = u[37]

    mod_level = clamp(sm / params.scale_sm, 0.0, 1.5)
    control_factor = if mode == direct_gmax
        1.0 + control_gain
    else
        1.0 + control_gain * mod_level
    end

    alphax_base = baseline_alphax(params, mode)
    g41_base = baseline_g41(params, mode)
    g32_base = baseline_g32(params, mode)

    alpha_exc_eff = mode == presynaptic ? alphax_base * control_factor : alphax_base
    g41_eff = (mode == postsynaptic || mode == direct_gmax) ? g41_base * control_factor : g41_base
    g32_eff = (mode == postsynaptic || mode == direct_gmax) ? g32_base * control_factor : g32_base

    Vs0 = 127.0 * V0 / 105.0 + 8265.0 / 105.0
    Vs1 = 127.0 * V1 / 105.0 + 8265.0 / 105.0
    Vs2 = 127.0 * V2 / 105.0 + 8265.0 / 105.0
    Vs3 = 127.0 * V3 / 105.0 + 8265.0 / 105.0
    Vs4 = 127.0 * V4 / 105.0 + 8265.0 / 105.0

    m0 = (0.1 * (50.0 - Vs0) / (exp((50.0 - Vs0) / 10.0) - 1.0)) /
         ((0.1 * (50.0 - Vs0) / (exp((50.0 - Vs0) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs0) / 18.0))
    m1 = (0.1 * (50.0 - Vs1) / (exp((50.0 - Vs1) / 10.0) - 1.0)) /
         ((0.1 * (50.0 - Vs1) / (exp((50.0 - Vs1) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs1) / 18.0))
    m2 = (0.1 * (50.0 - Vs2) / (exp((50.0 - Vs2) / 10.0) - 1.0)) /
         ((0.1 * (50.0 - Vs2) / (exp((50.0 - Vs2) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs2) / 18.0))
    m3 = (0.1 * (50.0 - Vs3) / (exp((50.0 - Vs3) / 10.0) - 1.0)) /
         ((0.1 * (50.0 - Vs3) / (exp((50.0 - Vs3) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs3) / 18.0))
    m4 = (0.1 * (50.0 - Vs4) / (exp((50.0 - Vs4) / 10.0) - 1.0)) /
         ((0.1 * (50.0 - Vs4) / (exp((50.0 - Vs4) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs4) / 18.0))

    applied_current = (p.protocol.current_on && t_ms >= params.t1_ms && t_ms <= params.t2_ms) ? params.Iapp : 0.0

    du[1] = applied_current +
            4.0 * m0^3 * h0 * (30.0 - V0) +
            0.3 * n0^4 * (-75.0 - V0) +
            0.01 * x0 * (30.0 - V0) +
            0.03 * Ca0 / (0.5 + Ca0) * (-75.0 - V0) +
            0.003 * (-40.0 - V0)

    du[2] = 4.0 * m1^3 * h1 * (30.0 - V1) +
            0.3 * n1^4 * (-75.0 - V1) +
            0.01 * x1 * (30.0 - V1) +
            0.03 * Ca1 / (0.5 + Ca1) * (-75.0 - V1) +
            0.003 * (-40.0 - V1) -
            g41_eff * (V1 - 30.0) * s4 -
            params.g21 * (V1 + 80.0) * s21 / params.scale2 +
            params.gelec * (V4 - V1)

    du[3] = 4.0 * m2^3 * h2 * (30.0 - V2) +
            0.3 * n2^4 * (-75.0 - V2) +
            0.01 * x2 * (30.0 - V2) +
            0.03 * Ca2 / (0.5 + Ca2) * (-75.0 - V2) +
            0.003 * (-40.0 - V2) -
            g32_eff * (V2 - 30.0) * s3 -
            params.g12 * (V2 + 80.0) * s12 / params.scale2 +
            params.gelec * (V3 - V2)

    du[4] = 4.0 * m3^3 * h3 * (30.0 - V3) +
            0.3 * n3^4 * (-75.0 - V3) +
            0.01 * x3 * (30.0 - V3) +
            0.03 * Ca3 / (0.5 + Ca3) * (-75.0 - V3) +
            0.003 * (-40.0 - V3) -
            params.g23 * (V3 + 80.0) * s2 / params.scalei -
            params.g43 * (V3 + 80.0) * s43 / params.scale3 +
            params.gelec * (V2 - V3) -
            params.g0 * (V3 - 30.0) * s0 / params.scale1

    du[5] = 4.0 * m4^3 * h4 * (30.0 - V4) +
            0.3 * n4^4 * (-75.0 - V4) +
            0.01 * x4 * (30.0 - V4) +
            0.03 * Ca4 / (0.5 + Ca4) * (-75.0 - V4) +
            0.003 * (-40.0 - V4) -
            params.g14 * (V4 + 80.0) * s1 / params.scalei -
            params.g34 * (V4 + 80.0) * s34 / params.scale3 +
            params.gelec * (V1 - V4) -
            params.g0 * (V4 - 30.0) * s0 / params.scale1

    du[11] = 0.0006 * (0.0085 * x0 * (140.0 - V0 + params.Ca_shift0) - Ca0)
    du[12] = 0.0003 * (0.0085 * x1 * (140.0 - V1 + params.Ca_shift1) - Ca1)
    du[13] = 0.0003 * (0.0085 * x2 * (140.0 - V2 + params.Ca_shift2) - Ca2)
    du[14] = 0.0003 * (0.0085 * x3 * (140.0 - V3 + params.Ca_shift3) - Ca3)
    du[15] = 0.0003 * (0.0085 * x4 * (140.0 - V4 + params.Ca_shift4) - Ca4)

    du[6] = ((1.0 / (exp(0.15 * (-V0 - 50.0 + params.x_shift)) + 1.0)) - x0) / 100.0
    du[7] = ((1.0 / (exp(0.15 * (-V1 - 50.0 + params.x_shift1)) + 1.0)) - x1) / 100.0
    du[8] = ((1.0 / (exp(0.15 * (-V2 - 50.0 + params.x_shift2)) + 1.0)) - x2) / 100.0
    du[9] = ((1.0 / (exp(0.15 * (-V3 - 50.0 + params.x_shift3)) + 1.0)) - x3) / 100.0
    du[10] = ((1.0 / (exp(0.15 * (-V4 - 50.0 + params.x_shift4)) + 1.0)) - x4) / 100.0

    du[16] = ((1.0 - h0) * (0.07 * exp((25.0 - Vs0) / 20.0)) - h0 * (1.0 / (1.0 + exp((55.0 - Vs0) / 10.0)))) / 12.5
    du[17] = ((1.0 - h1) * (0.07 * exp((25.0 - Vs1) / 20.0)) - h1 * (1.0 / (1.0 + exp((55.0 - Vs1) / 10.0)))) / 12.5
    du[18] = ((1.0 - h2) * (0.07 * exp((25.0 - Vs2) / 20.0)) - h2 * (1.0 / (1.0 + exp((55.0 - Vs2) / 10.0)))) / 12.5
    du[19] = ((1.0 - h3) * (0.07 * exp((25.0 - Vs3) / 20.0)) - h3 * (1.0 / (1.0 + exp((55.0 - Vs3) / 10.0)))) / 12.5
    du[20] = ((1.0 - h4) * (0.07 * exp((25.0 - Vs4) / 20.0)) - h4 * (1.0 / (1.0 + exp((55.0 - Vs4) / 10.0)))) / 12.5

    du[21] = ((1.0 - n0) * (0.01 * (55.0 - Vs0) / (exp((55.0 - Vs0) / 10.0) - 1.0)) - n0 * (0.125 * exp((45.0 - Vs0) / 80.0))) / 12.5
    du[22] = ((1.0 - n1) * (0.01 * (55.0 - Vs1) / (exp((55.0 - Vs1) / 10.0) - 1.0)) - n1 * (0.125 * exp((45.0 - Vs1) / 80.0))) / 12.5
    du[23] = ((1.0 - n2) * (0.01 * (55.0 - Vs2) / (exp((55.0 - Vs2) / 10.0) - 1.0)) - n2 * (0.125 * exp((45.0 - Vs2) / 80.0))) / 12.5
    du[24] = ((1.0 - n3) * (0.01 * (55.0 - Vs3) / (exp((55.0 - Vs3) / 10.0) - 1.0)) - n3 * (0.125 * exp((45.0 - Vs3) / 80.0))) / 12.5
    du[25] = ((1.0 - n4) * (0.01 * (55.0 - Vs4) / (exp((55.0 - Vs4) / 10.0) - 1.0)) - n4 * (0.125 * exp((45.0 - Vs4) / 80.0))) / 12.5

    du[26] = 0.5 * ((1.0 / (1.0 + exp(10.0 * (V1 + 53.0)))) - y1) / (7.1 + 10.4 / (1.0 + exp((V1 + 68.0) / 2.2)))
    du[27] = 0.5 * ((1.0 / (1.0 + exp(10.0 * (V2 + 53.0)))) - y2) / (7.1 + 10.4 / (1.0 + exp((V2 + 68.0) / 2.2)))

    du[28] = params.alpha1 * (1.0 - s0) * syn_activation(V0) - params.beta1 * (s0 - params.s0_floor)
    du[29] = params.alphai * s1 * (1.0 - s1) * syn_activation(V1) - params.betai * (s1 - params.si_inhib_floor)
    du[30] = params.alphai * s2 * (1.0 - s2) * syn_activation(V2) - params.betai * (s2 - params.si_inhib_floor)
    du[31] = alpha_exc_eff * s3 * (1.0 - s3) * syn_activation(V3) - params.betax * (s3 - params.si3_exc_floor)
    du[32] = alpha_exc_eff * s4 * (1.0 - s4) * syn_activation(V4) - params.betax * (s4 - params.si3_exc_floor)

    du[33] = params.alpha2 * (1.0 - s12) * syn_activation(V1) - params.beta2 * s12
    du[34] = params.alpha2 * (1.0 - s21) * syn_activation(V2) - params.beta2 * s21
    du[35] = params.alpha3 * (1.0 - s34) * syn_activation(V3) - params.beta3 * s34
    du[36] = params.alpha3 * (1.0 - s43) * syn_activation(V4) - params.beta3 * s43

    du[37] = params.alpham * sm * (1.0 - sm) * syn_activation(V0) - params.betam * (sm - params.sm_floor)
    return nothing
end

function run_simulation(params::ModelParams, control_gain::Float64, mode::ControlMode, config::SweepConfig; terminate_on_withdrawal_end::Bool = true, initial_u0 = nothing)
    sim_params = mode_params(params, mode)
    required_complete_bursts = config.transient_bursts + config.measured_cycles
    tracker = SpikeTracker(Float64[], -Inf, config.spike_refractory_s * 1000.0, config.spike_threshold_mv, required_complete_bursts)
    protocol = ProtocolTracker(:seeking_onset, true, false, NaN, NaN, NaN, NaN)
    ode_params = (params = sim_params, control_gain = control_gain, mode = mode, config = config, protocol = protocol)
    spike_callback = make_spike_callback(tracker)
    protocol_callback = make_protocol_callback(protocol, terminate_on_withdrawal_end)
    function spike_affect!(integrator)
        t_ms = integrator.t
        if isempty(tracker.spike_times_ms) || (t_ms - tracker.last_spike_ms) >= tracker.spike_refractory_ms
            push!(tracker.spike_times_ms, t_ms)
            tracker.last_spike_ms = t_ms
            update_protocol_from_spikes!(protocol, tracker.spike_times_ms, config, integrator)
        end
    end
    spike_callback = ContinuousCallback(
        (u, t, integrator) -> u[2] - tracker.threshold_mv,
        spike_affect!;
        affect_neg! = nothing,
        save_positions = (false, false),
    )
    callback = CallbackSet(spike_callback, protocol_callback)
    tspan_ms = (0.0, config.max_time_s * 1000.0)
    u0 = isnothing(initial_u0) ? initial_state() : copy(initial_u0)
    problem = ODEProblem(network_ode!, u0, tspan_ms, ode_params)
    solution = solve(
        problem,
        RK4();
        saveat = config.saveat_ms,
        reltol = 1e-6,
        abstol = 1e-6,
        callback = callback,
        maxiters = 10^9,
    )

    time_s = solution.t ./ 1000.0
    states = Array(solution)

    return (
        time_s = time_s,
        V0 = states[1, :],
        V1 = states[2, :],
        V2 = states[3, :],
        V3 = states[4, :],
        V4 = states[5, :],
        states = states,
        s3 = states[31, :],
        s4 = states[32, :],
        sm = states[37, :],
        spike_times_s = tracker.spike_times_ms ./ 1000.0,
        params = sim_params,
        onset_failed = protocol.onset_failed,
        stimulus_off_s = protocol.stimulus_off_s,
        withdrawal_end_s = protocol.withdrawal_end_s,
        measured_start_s = protocol.measured_start_s,
        measured_stop_s = protocol.measured_stop_s,
        current_on = protocol.current_on,
    )
end

function effective_excitation_trace(sm, s_exc, control_gain, mode::ControlMode, params::ModelParams)
    mod_level = clamp.(sm ./ params.scale_sm, 0.0, 1.5)
    gain = if mode == direct_gmax
        fill(1.0 + control_gain, length(mod_level))
    else
        1.0 .+ control_gain .* mod_level
    end

    g_base = baseline_g41(params, mode)
    alpha_base = baseline_alphax(params, mode)
    g_eff = (mode == postsynaptic || mode == direct_gmax) ? g_base .* gain : fill(g_base, length(gain))
    alpha_eff = mode == presynaptic ? alpha_base .* gain : fill(alpha_base, length(gain))
    return (
        strength = g_eff .* s_exc,
        gain = gain,
        alpha_eff = alpha_eff,
        g_eff = g_eff,
    )
end

function interval_indices(time_s::AbstractVector{<:Real}, t_start::Float64, t_stop::Float64)
    idx = findall(t -> t >= t_start && t <= t_stop, time_s)
    if !isempty(idx)
        return idx
    end

    start_idx = clamp(searchsortedfirst(time_s, t_start), 1, length(time_s))
    stop_idx = clamp(searchsortedlast(time_s, t_stop), 1, length(time_s))
    if stop_idx < start_idx
        stop_idx = start_idx
    end
    return start_idx:stop_idx
end

function collect_burst_metrics(sim, control_gain::Float64, mode::ControlMode, params::ModelParams, config::SweepConfig)
    required_bursts = config.transient_bursts + config.measured_cycles
    bursts = detect_bursts_from_spikes(sim.spike_times_s, config)
    complete_start_times = complete_burst_start_times(bursts)
    if sim.onset_failed || length(complete_start_times) < required_bursts
        return DataFrame(
            mode = String[],
            control_gain = Float64[],
            burst_time_s = Float64[],
            burst_frequency_hz = Float64[],
            mean_effective_excitation = Float64[],
            peak_effective_excitation = Float64[],
            mean_gain_factor = Float64[],
            alpha_eff = Float64[],
            g_eff = Float64[],
        )
    end

    traces = effective_excitation_trace(sim.sm, sim.s4, control_gain, mode, sim.params)
    rows = DataFrame(
        mode = String[],
        control_gain = Float64[],
        burst_time_s = Float64[],
        burst_frequency_hz = Float64[],
        mean_effective_excitation = Float64[],
        peak_effective_excitation = Float64[],
        mean_gain_factor = Float64[],
        alpha_eff = Float64[],
        g_eff = Float64[],
    )

    measured_start_times = complete_start_times[(config.transient_bursts + 1):(config.transient_bursts + config.measured_cycles)]
    for k in 1:(length(measured_start_times) - 1)
        interval_start = measured_start_times[k]
        interval_stop = measured_start_times[k + 1]
        burst_time = interval_start
        ibp = interval_stop - interval_start
        if ibp <= 0
            continue
        end

        interval_idx = interval_indices(sim.time_s, interval_start, interval_stop)
        push!(rows, (
            mode_label(mode),
            control_gain,
            burst_time,
            1.0 / ibp,
            mean(traces.strength[interval_idx]),
            maximum(traces.strength[interval_idx]),
            mean(traces.gain[interval_idx]),
            mean(traces.alpha_eff[interval_idx]),
            mean(traces.g_eff[interval_idx]),
        ))
    end

    return rows
end

function post_withdrawal_metrics(sim, config::SweepConfig)
    if sim.onset_failed || isnan(sim.stimulus_off_s) || isnan(sim.withdrawal_end_s)
        return (count = NaN, mean_frequency_hz = NaN)
    end

    bursts = detect_bursts_from_spikes(sim.spike_times_s, config)
    complete_starts = complete_burst_start_times(bursts)
    complete_ends = complete_burst_end_times(bursts)
    keep = findall(i -> complete_starts[i] >= sim.stimulus_off_s && complete_ends[i] <= sim.withdrawal_end_s, eachindex(complete_starts))
    post_starts = complete_starts[keep]
    post_count = length(post_starts)
    post_mean_frequency = post_count >= 2 ? mean(1.0 ./ diff(post_starts)) : NaN
    return (count = Float64(post_count), mean_frequency_hz = post_mean_frequency)
end

function summarize_results(points::DataFrame)
    if isempty(points)
        return DataFrame(
            mode = String[],
            control_gain = Float64[],
            onset_failed = Bool[],
            n_pre_intervals = Int[],
            mean_pre_burst_frequency_hz = Float64[],
            std_pre_burst_frequency_hz = Float64[],
            post_withdrawal_burst_count = Float64[],
            post_withdrawal_mean_frequency_hz = Float64[],
            post_pre_frequency_ratio = Float64[],
            stimulus_off_s = Float64[],
            withdrawal_duration_s = Float64[],
            mean_alpha_eff = Float64[],
            mean_g_eff = Float64[],
        )
    end

    grouped = groupby(points, [:mode, :control_gain])
    return combine(grouped,
        nrow => :n_pre_intervals,
        :burst_frequency_hz => mean => :mean_pre_burst_frequency_hz,
        :burst_frequency_hz => std => :std_pre_burst_frequency_hz,
        :alpha_eff => mean => :mean_alpha_eff,
        :g_eff => mean => :mean_g_eff,
    )
end

function summarize_single_run(points::DataFrame, sim, control_gain::Float64, mode::ControlMode, params::ModelParams, config::SweepConfig)
    post = post_withdrawal_metrics(sim, config)
    withdrawal_duration_s = isnan(sim.withdrawal_end_s) || isnan(sim.stimulus_off_s) ? NaN : sim.withdrawal_end_s - sim.stimulus_off_s

    if isempty(points)
        return DataFrame(
            mode = [mode_label(mode)],
            control_gain = [control_gain],
            onset_failed = [sim.onset_failed],
            n_pre_intervals = [0],
            mean_pre_burst_frequency_hz = [NaN],
            std_pre_burst_frequency_hz = [NaN],
            post_withdrawal_burst_count = [post.count],
            post_withdrawal_mean_frequency_hz = [post.mean_frequency_hz],
            post_pre_frequency_ratio = [NaN],
            stimulus_off_s = [sim.stimulus_off_s],
            withdrawal_duration_s = [withdrawal_duration_s],
            mean_alpha_eff = [NaN],
            mean_g_eff = [NaN],
        )
    end

    pre_mean_frequency = mean(points.burst_frequency_hz)
    post_pre_ratio = isfinite(post.mean_frequency_hz) && isfinite(pre_mean_frequency) && pre_mean_frequency > 0 ? post.mean_frequency_hz / pre_mean_frequency : NaN

    return DataFrame(
        mode = [mode_label(mode)],
        control_gain = [control_gain],
        onset_failed = [sim.onset_failed],
        n_pre_intervals = [nrow(points)],
        mean_pre_burst_frequency_hz = [pre_mean_frequency],
        std_pre_burst_frequency_hz = [coalesce(std(points.burst_frequency_hz), 0.0)],
        post_withdrawal_burst_count = [post.count],
        post_withdrawal_mean_frequency_hz = [post.mean_frequency_hz],
        post_pre_frequency_ratio = [post_pre_ratio],
        stimulus_off_s = [sim.stimulus_off_s],
        withdrawal_duration_s = [withdrawal_duration_s],
        mean_alpha_eff = [mean(points.alpha_eff)],
        mean_g_eff = [mean(points.g_eff)],
    )
end

function validate_summary(summary::DataFrame)
    for mode_name in ("direct_gmax", "presynaptic", "postsynaptic")
        subset = sort(summary[summary.mode .== mode_name, :], :control_gain)
        isempty(subset) && continue
        subset.control_gain[1] == 0.0 || error("Mode $(mode_name) sweep no longer begins at zero gain")
    end
end

function spec_lookup(specs::Vector{ModeSweepSpec})
    Dict(mode_label(spec.mode) => spec for spec in specs)
end

function plot_results(points::DataFrame, summary::DataFrame, specs::Vector{ModeSweepSpec})
    colors = Dict(
        "direct_gmax" => :darkgreen,
        "presynaptic" => :dodgerblue3,
        "postsynaptic" => :firebrick3,
    )
    ordered_modes = ("direct_gmax", "presynaptic", "postsynaptic")

    p1 = plot(
        title = "Burst Frequency vs Control Gain",
        xlabel = "Peak Si3→Si2 Excitation Per Cycle (g*s proxy)",
        ylabel = "Mean Burst Frequency (Hz)",
        legend = :topleft,
        framestyle = :box,
        grid = true,
    )
    xlabel!(p1, "Control Gain")

    for mode_name in ordered_modes
        subset = summary[summary.mode .== mode_name, :]
        isempty(subset) && continue
        subset = sort(subset, :control_gain)
        ribbon = coalesce.(subset.std_burst_frequency_hz, 0.0)
        plot!(
            p1,
            subset.control_gain,
            subset.mean_burst_frequency_hz;
            ribbon = ribbon,
            label = mode_name,
            color = colors[mode_name],
            linewidth = 2.5,
            marker = :circle,
            markersize = 4,
            fillalpha = 0.12,
        )
    end

    p2 = plot(
        title = "Effective Excitation vs Control Gain",
        xlabel = "Control Gain",
        ylabel = "Mean Peak Si3->Si2 Excitation Per Cycle",
        legend = :topleft,
        framestyle = :box,
        grid = true,
    )

    for mode_name in ordered_modes
        subset = summary[summary.mode .== mode_name, :]
        isempty(subset) && continue
        subset = sort(subset, :control_gain)
        plot!(
            p2,
            subset.control_gain,
            subset.peak_cycle_excitation;
            label = mode_name,
            color = colors[mode_name],
            linewidth = 2.5,
            marker = :circle,
            markersize = 4,
        )
    end

    p3 = scatter(
        title = "Burst Frequency vs Effective Excitation",
        xlabel = "Peak Si3->Si2 Excitation Per Cycle (g*s proxy)",
        ylabel = "Burst Frequency (Hz)",
        legend = :topleft,
        framestyle = :box,
        grid = true,
    )

    for mode_name in ordered_modes
        subset = points[points.mode .== mode_name, :]
        isempty(subset) && continue
        scatter!(
            p3,
            subset.peak_effective_excitation,
            subset.burst_frequency_hz;
            label = mode_name,
            color = colors[mode_name],
            alpha = 0.6,
            markersize = 5,
            markerstrokewidth = 0.3,
        )
    end

    return plot(
        p1,
        p2,
        p3;
        layout = (1, 3),
        size = (1900, 600),
        bottom_margin = 10mm,
        left_margin = 4mm,
        right_margin = 4mm,
        top_margin = 4mm,
    )
end

function plot_mode_panels(raw_points::DataFrame, summary::DataFrame, specs::Vector{ModeSweepSpec})
    colors = Dict(
        "direct_gmax" => :darkgreen,
        "presynaptic" => :dodgerblue3,
        "postsynaptic" => :firebrick3,
    )
    panels = Plots.Plot[]

    for spec in specs
        mode_name = mode_label(spec.mode)
        points_subset = raw_points[raw_points.mode .== mode_name, :]
        summary_subset = sort(summary[summary.mode .== mode_name, :], :control_gain)
        p_freq_control = plot(
            title = "$(spec.label): Pre Frequency vs Control",
            xlabel = spec.control_axis_label,
            ylabel = "Mean Pre-Withdrawal Burst Frequency (Hz)",
            legend = false,
            framestyle = :box,
            grid = true,
        )
        if !isempty(points_subset)
            scatter!(
                p_freq_control,
                points_subset.control_gain,
                points_subset.burst_frequency_hz;
                color = colors[mode_name],
                alpha = 0.04,
                markersize = 2.5,
                markerstrokewidth = 0.1,
            )
        end
        if !isempty(summary_subset)
            scatter!(
                p_freq_control,
                summary_subset.control_gain,
                summary_subset.mean_pre_burst_frequency_hz;
                color = colors[mode_name],
                alpha = 0.95,
                markersize = 7,
                markerstrokewidth = 0.4,
            )
        end

        p_post_count = plot(
            title = "$(spec.label): Post-Withdrawal Bursts vs Control",
            xlabel = spec.control_axis_label,
            ylabel = "Post-Withdrawal Complete Burst Count",
            legend = false,
            framestyle = :box,
            grid = true,
        )
        if !isempty(summary_subset)
            scatter!(
                p_post_count,
                summary_subset.control_gain,
                summary_subset.post_withdrawal_burst_count;
                color = colors[mode_name],
                alpha = 0.95,
                markersize = 7,
                markerstrokewidth = 0.4,
            )
        end

        p_post_ratio = plot(
            title = "$(spec.label): Post/Pre Ratio vs Control",
            xlabel = spec.control_axis_label,
            ylabel = "Post/Pre Burst Frequency Ratio",
            legend = false,
            framestyle = :box,
            grid = true,
        )
        if !isempty(summary_subset)
            scatter!(
                p_post_ratio,
                summary_subset.control_gain,
                summary_subset.post_pre_frequency_ratio;
                color = colors[mode_name],
                alpha = 0.95,
                markersize = 7,
                markerstrokewidth = 0.4,
            )
        end

        append!(panels, [p_freq_control, p_post_count, p_post_ratio])
    end

    return plot(
        panels...;
        layout = (length(specs), 3),
        size = (1800, 1500),
        bottom_margin = 10mm,
        left_margin = 6mm,
        right_margin = 4mm,
        top_margin = 6mm,
    )
end

function generate_data(config::SweepConfig, params::ModelParams, specs::Vector{ModeSweepSpec})
    raw_points = DataFrame(
        mode = String[],
        control_gain = Float64[],
        burst_time_s = Float64[],
        burst_frequency_hz = Float64[],
        mean_effective_excitation = Float64[],
        peak_effective_excitation = Float64[],
        mean_gain_factor = Float64[],
        alpha_eff = Float64[],
        g_eff = Float64[],
    )
    summary = DataFrame(
        mode = String[],
        control_gain = Float64[],
        onset_failed = Bool[],
        n_pre_intervals = Int[],
        mean_pre_burst_frequency_hz = Float64[],
        std_pre_burst_frequency_hz = Float64[],
        post_withdrawal_burst_count = Float64[],
        post_withdrawal_mean_frequency_hz = Float64[],
        post_pre_frequency_ratio = Float64[],
        stimulus_off_s = Float64[],
        withdrawal_duration_s = Float64[],
        mean_alpha_eff = Float64[],
        mean_g_eff = Float64[],
    )

    for spec in specs
        println("Running $(spec.label) sweep...")
        for control_gain in spec.control_values
            println("  control_gain = $(control_gain)")
            sim = run_simulation(params, control_gain, spec.mode, config)
            points = collect_burst_metrics(sim, control_gain, spec.mode, params, config)
            append!(raw_points, points)
            append!(summary, summarize_single_run(points, sim, control_gain, spec.mode, params, config))
            println("    collected $(nrow(points)) burst intervals")
        end
    end

    return raw_points, summary
end

function save_data(raw_points::DataFrame, summary::DataFrame)
    CSV.write(RAW_CSV, raw_points)
    CSV.write(SUMMARY_CSV, summary)
end

function load_summary()
    isfile(SUMMARY_CSV) || error("Summary CSV not found at $(SUMMARY_CSV). Run with `generate` first.")
    return CSV.read(SUMMARY_CSV, DataFrame)
end

function load_raw_points()
    isfile(RAW_CSV) || error("Raw points CSV not found at $(RAW_CSV). Run with `generate` first.")
    return CSV.read(RAW_CSV, DataFrame)
end

function save_figure(raw_points::DataFrame, summary::DataFrame, specs::Vector{ModeSweepSpec})
    validate_summary(summary)
    figure = plot_mode_panels(raw_points, summary, specs)
    savefig(figure, FIGURE_PNG)
end

function main(args = ARGS)
    mkpath(OUTPUT_DIR)

    config = default_config()
    params = default_params()
    specs = default_sweep_specs()
    mode = isempty(args) ? "both" : lowercase(args[1])

    if mode == "generate"
        raw_points, summary = generate_data(config, params, specs)
        save_data(raw_points, summary)
        println()
        println("Saved data:")
        println("  $(RAW_CSV)")
        println("  $(SUMMARY_CSV)")
        return
    end

    if mode == "plot"
        raw_points = load_raw_points()
        summary = load_summary()
        save_figure(raw_points, summary, specs)
        println()
        println("Saved figure:")
        println("  $(FIGURE_PNG)")
        return
    end

    if mode == "both"
        raw_points, summary = generate_data(config, params, specs)
        save_data(raw_points, summary)
        save_figure(raw_points, summary, specs)
        println()
        println("Saved:")
        println("  $(RAW_CSV)")
        println("  $(SUMMARY_CSV)")
        println("  $(FIGURE_PNG)")
        return
    end

    error("Unknown mode `$(mode)`. Use `generate`, `plot`, or `both`.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
