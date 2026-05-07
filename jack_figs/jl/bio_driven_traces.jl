using CSV
using DataFrames
using OrdinaryDiffEq

include(joinpath(@__DIR__, "calibrated_neuromod.jl"))

const BIO_OUTPUT_DIR = joinpath(@__DIR__, "outputs")
const AKIRA_DIR = joinpath(dirname(@__DIR__), "akira", "extracted")
const BIO_SI1_TXT = joinpath(AKIRA_DIR, "130618-03 Si1.txt")
const BIO_SI2_TXT = joinpath(AKIRA_DIR, "130618-03 Si2.txt")
const BIO_DRIVEN_TRACES_CSV = joinpath(BIO_OUTPUT_DIR, "bio_driven_traces.csv")
const BIO_DRIVEN_POINTS_CSV = joinpath(BIO_OUTPUT_DIR, "bio_driven_points.csv")

const PRESYNAPTIC_TRACE_PARAMS = calibrated_params(presynaptic)
const POSTSYNAPTIC_TRACE_PARAMS = calibrated_params(postsynaptic)
const PRESYNAPTIC_TRACE_GAIN = PRESYNAPTIC_ANCHOR_GAIN
const POSTSYNAPTIC_TRACE_GAIN = POSTSYNAPTIC_ANCHOR_GAIN
const TRACE_SAVEAT_MS = 10.0

function load_two_column_txt(path::AbstractString)
    df = DataFrame(CSV.File(path; delim = '\t', normalizenames = false))
    return Float64.(df[!, 1]), Float64.(df[!, 2])
end

function interpolate_drive(time_s::Float64, drive_t::Vector{Float64}, drive_v::Vector{Float64})
    time_s <= drive_t[1] && return drive_v[1]
    time_s >= drive_t[end] && return drive_v[end]
    idx = clamp(searchsortedlast(drive_t, time_s), 1, length(drive_t) - 1)
    t0 = drive_t[idx]
    t1 = drive_t[idx + 1]
    v0 = drive_v[idx]
    v1 = drive_v[idx + 1]
    frac = (time_s - t0) / max(t1 - t0, eps())
    return v0 + frac * (v1 - v0)
end

function publication_driven_network_ode!(du, u, p, t_ms)
    params = p.params
    control_gain = p.control_gain
    mode = p.mode
    V0 = interpolate_drive(t_ms / 1000.0, p.drive_t, p.drive_v)

    _, V1, V2, V3, V4 = u[1:5]
    _, x1, x2, x3, x4 = u[6:10]
    _, Ca1, Ca2, Ca3, Ca4 = u[11:15]
    _, h1, h2, h3, h4 = u[16:20]
    _, n1, n2, n3, n4 = u[21:25]
    y1, y2 = u[26:27]
    s0, s1, s2, s3, s4 = u[28:32]
    s12, s21, s34, s43 = u[33:36]
    sm = u[37]

    mod_level = clamp(sm / params.scale_sm, 0.0, 1.5)
    control_factor = 1.0 + control_gain * mod_level

    alphax_base = baseline_alphax(params, mode)
    g41_base = baseline_g41(params, mode)
    g32_base = baseline_g32(params, mode)
    alpha_exc_eff = mode == presynaptic ? alphax_base * control_factor : alphax_base
    g41_eff = mode == postsynaptic ? g41_base * control_factor : g41_base
    g32_eff = mode == postsynaptic ? g32_base * control_factor : g32_base

    Vs1 = 127.0 * V1 / 105.0 + 8265.0 / 105.0
    Vs2 = 127.0 * V2 / 105.0 + 8265.0 / 105.0
    Vs3 = 127.0 * V3 / 105.0 + 8265.0 / 105.0
    Vs4 = 127.0 * V4 / 105.0 + 8265.0 / 105.0

    m1 = (0.1 * (50.0 - Vs1) / (exp((50.0 - Vs1) / 10.0) - 1.0)) / ((0.1 * (50.0 - Vs1) / (exp((50.0 - Vs1) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs1) / 18.0))
    m2 = (0.1 * (50.0 - Vs2) / (exp((50.0 - Vs2) / 10.0) - 1.0)) / ((0.1 * (50.0 - Vs2) / (exp((50.0 - Vs2) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs2) / 18.0))
    m3 = (0.1 * (50.0 - Vs3) / (exp((50.0 - Vs3) / 10.0) - 1.0)) / ((0.1 * (50.0 - Vs3) / (exp((50.0 - Vs3) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs3) / 18.0))
    m4 = (0.1 * (50.0 - Vs4) / (exp((50.0 - Vs4) / 10.0) - 1.0)) / ((0.1 * (50.0 - Vs4) / (exp((50.0 - Vs4) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs4) / 18.0))

    du[1] = 0.0
    du[2] = 4.0 * m1^3 * h1 * (30.0 - V1) + 0.3 * n1^4 * (-75.0 - V1) + 0.01 * x1 * (30.0 - V1) + 0.03 * Ca1 / (0.5 + Ca1) * (-75.0 - V1) + 0.003 * (-40.0 - V1) - g41_eff * (V1 - 30.0) * s4 - params.g21 * (V1 + 80.0) * s21 / params.scale2 + params.gelec * (V4 - V1)
    du[3] = 4.0 * m2^3 * h2 * (30.0 - V2) + 0.3 * n2^4 * (-75.0 - V2) + 0.01 * x2 * (30.0 - V2) + 0.03 * Ca2 / (0.5 + Ca2) * (-75.0 - V2) + 0.003 * (-40.0 - V2) - g32_eff * (V2 - 30.0) * s3 - params.g12 * (V2 + 80.0) * s12 / params.scale2 + params.gelec * (V3 - V2)
    du[4] = 4.0 * m3^3 * h3 * (30.0 - V3) + 0.3 * n3^4 * (-75.0 - V3) + 0.01 * x3 * (30.0 - V3) + 0.03 * Ca3 / (0.5 + Ca3) * (-75.0 - V3) + 0.003 * (-40.0 - V3) - params.g23 * (V3 + 80.0) * s2 / params.scalei - params.g43 * (V3 + 80.0) * s43 / params.scale3 + params.gelec * (V2 - V3) - params.g0 * (V3 - 30.0) * s0 / params.scale1
    du[5] = 4.0 * m4^3 * h4 * (30.0 - V4) + 0.3 * n4^4 * (-75.0 - V4) + 0.01 * x4 * (30.0 - V4) + 0.03 * Ca4 / (0.5 + Ca4) * (-75.0 - V4) + 0.003 * (-40.0 - V4) - params.g14 * (V4 + 80.0) * s1 / params.scalei - params.g34 * (V4 + 80.0) * s34 / params.scale3 + params.gelec * (V1 - V4) - params.g0 * (V4 - 30.0) * s0 / params.scale1

    du[6] = 0.0
    du[7] = ((1.0 / (exp(0.15 * (-V1 - 50.0 + params.x_shift1)) + 1.0)) - x1) / 100.0
    du[8] = ((1.0 / (exp(0.15 * (-V2 - 50.0 + params.x_shift2)) + 1.0)) - x2) / 100.0
    du[9] = ((1.0 / (exp(0.15 * (-V3 - 50.0 + params.x_shift3)) + 1.0)) - x3) / 100.0
    du[10] = ((1.0 / (exp(0.15 * (-V4 - 50.0 + params.x_shift4)) + 1.0)) - x4) / 100.0

    du[11] = 0.0
    du[12] = 0.0003 * (0.0085 * x1 * (140.0 - V1 + params.Ca_shift1) - Ca1)
    du[13] = 0.0003 * (0.0085 * x2 * (140.0 - V2 + params.Ca_shift2) - Ca2)
    du[14] = 0.0003 * (0.0085 * x3 * (140.0 - V3 + params.Ca_shift3) - Ca3)
    du[15] = 0.0003 * (0.0085 * x4 * (140.0 - V4 + params.Ca_shift4) - Ca4)

    du[16] = 0.0
    du[17] = ((1.0 - h1) * (0.07 * exp((25.0 - Vs1) / 20.0)) - h1 * (1.0 / (1.0 + exp((55.0 - Vs1) / 10.0)))) / 12.5
    du[18] = ((1.0 - h2) * (0.07 * exp((25.0 - Vs2) / 20.0)) - h2 * (1.0 / (1.0 + exp((55.0 - Vs2) / 10.0)))) / 12.5
    du[19] = ((1.0 - h3) * (0.07 * exp((25.0 - Vs3) / 20.0)) - h3 * (1.0 / (1.0 + exp((55.0 - Vs3) / 10.0)))) / 12.5
    du[20] = ((1.0 - h4) * (0.07 * exp((25.0 - Vs4) / 20.0)) - h4 * (1.0 / (1.0 + exp((55.0 - Vs4) / 10.0)))) / 12.5

    du[21] = 0.0
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

function representative_trace_initial_state(params::ModelParams, gain::Float64, mode::ControlMode)
    return settled_initial_state(params, gain, mode)
end

function simulate_bio_driven_trace(params::ModelParams, gain::Float64, mode::ControlMode, drive_time_s::Vector{Float64}, drive_voltage_mv::Vector{Float64})
    sim_params = mode_params(params, mode)
    u0 = representative_trace_initial_state(params, gain, mode)
    problem = ODEProblem(
        publication_driven_network_ode!,
        u0,
        (0.0, last(drive_time_s) * 1000.0),
        (params = sim_params, control_gain = gain, mode = mode, drive_t = drive_time_s, drive_v = drive_voltage_mv),
    )
    solution = solve(problem, RK4(); saveat = TRACE_SAVEAT_MS, reltol = 1e-6, abstol = 1e-6, maxiters = 10^9)
    states = Array(solution)
    time_s = solution.t ./ 1000.0
    return DataFrame(
        mode = fill(mode_label(mode), length(time_s)),
        time_s = time_s,
        V0 = [interpolate_drive(t, drive_time_s, drive_voltage_mv) for t in time_s],
        V1 = states[2, :],
    )
end

function threshold_crossing_spikes(time_s::AbstractVector{<:Real}, voltage_mv::AbstractVector{<:Real}; threshold_mv::Float64 = 0.0, refractory_s::Float64 = 0.02)
    spikes = Float64[]
    last_spike = -Inf
    n = min(length(time_s), length(voltage_mv))
    for i in 2:n
        v_prev = Float64(voltage_mv[i - 1])
        v_now = Float64(voltage_mv[i])
        t_prev = Float64(time_s[i - 1])
        t_now = Float64(time_s[i])
        if v_prev < threshold_mv && v_now >= threshold_mv
            frac = (threshold_mv - v_prev) / max(v_now - v_prev, eps())
            spike_t = t_prev + frac * (t_now - t_prev)
            if isempty(spikes) || (spike_t - last_spike) >= refractory_s
                push!(spikes, spike_t)
                last_spike = spike_t
            end
        end
    end
    return spikes
end

function spike_burst_config(; spike_threshold_mv::Float64 = -20.0, spike_refractory_s::Float64 = 0.02, burst_factor::Float64 = 2.0)
    base = calibrated_scan_config()
    return SweepConfig(
        base.max_time_s,
        base.saveat_ms,
        spike_threshold_mv,
        spike_refractory_s,
        burst_factor,
        base.onset_timeout_s,
        base.transient_bursts,
        base.measured_cycles,
    )
end

function nearest_time_indices(time_s::AbstractVector{<:Real}, event_times_s::AbstractVector{<:Real})
    indices = Int[]
    isempty(time_s) && return indices
    for t in Float64.(event_times_s)
        idx = clamp(searchsortedfirst(time_s, t), 1, length(time_s))
        if idx > 1 && abs(Float64(time_s[idx - 1]) - t) < abs(Float64(time_s[idx]) - t)
            idx -= 1
        end
        push!(indices, idx)
    end
    return indices
end

function detect_bursts_from_trace_spikes(time_s::AbstractVector{<:Real}, voltage_mv::AbstractVector{<:Real}; spike_threshold_mv::Float64 = -20.0, spike_refractory_s::Float64 = 0.02, burst_factor::Float64 = 2.0)
    spike_times = threshold_crossing_spikes(time_s, voltage_mv; threshold_mv = spike_threshold_mv, refractory_s = spike_refractory_s)
    bursts = detect_bursts_from_spikes(spike_times, spike_burst_config(; spike_threshold_mv = spike_threshold_mv, spike_refractory_s = spike_refractory_s, burst_factor = burst_factor))
    starts = nearest_time_indices(time_s, complete_burst_start_times(bursts))
    ends = nearest_time_indices(time_s, complete_burst_end_times(bursts))
    return starts, ends
end

function matched_si1_vs_si2_points(trace::DataFrame; spike_threshold_mv::Float64 = -20.0, spike_refractory_s::Float64 = 0.02, burst_factor::Float64 = 2.0, skip_initial_isis::Int = 3)
    time_s = trace[!, :time_s]
    v0 = trace[!, :V0]
    v1 = trace[!, :V1]
    burst_starts, _ = detect_bursts_from_trace_spikes(time_s, v1; spike_threshold_mv = spike_threshold_mv, spike_refractory_s = spike_refractory_s, burst_factor = burst_factor)
    burst_times = time_s[burst_starts]
    if length(burst_times) < 2
        return DataFrame(mode = String[], si1_frequency_hz = Float64[], burst_frequency_hz = Float64[])
    end

    spike_times = threshold_crossing_spikes(time_s, v0; threshold_mv = spike_threshold_mv, refractory_s = spike_refractory_s)
    if length(spike_times) < 2
        return DataFrame(mode = String[], si1_frequency_hz = Float64[], burst_frequency_hz = Float64[])
    end

    burst_freq = 1.0 ./ diff(burst_times)

    si1_at_bursts = Float64[]
    for (t_start, t_stop) in zip(burst_times[1:end-1], burst_times[2:end])
        n_spikes = count(t -> t_start <= t < t_stop, spike_times)
        push!(si1_at_bursts, n_spikes / max(t_stop - t_start, eps()))
    end
    keep = (max(skip_initial_isis, 0) + 1):length(burst_freq)
    isempty(keep) && return DataFrame(mode = String[], si1_frequency_hz = Float64[], burst_frequency_hz = Float64[])

    return DataFrame(
        mode = fill(String(trace[1, :mode]), length(keep)),
        si1_frequency_hz = si1_at_bursts[keep],
        burst_frequency_hz = burst_freq[keep],
    )
end

function main()
    bio_t1, bio_v1 = load_two_column_txt(BIO_SI1_TXT)
    bio_t1 .-= first(bio_t1)
    traces = vcat(
        simulate_bio_driven_trace(PRESYNAPTIC_TRACE_PARAMS, PRESYNAPTIC_TRACE_GAIN, presynaptic, bio_t1, bio_v1),
        simulate_bio_driven_trace(POSTSYNAPTIC_TRACE_PARAMS, POSTSYNAPTIC_TRACE_GAIN, postsynaptic, bio_t1, bio_v1),
    )
    points = vcat(
        matched_si1_vs_si2_points(traces[traces.mode .== "presynaptic", :]),
        matched_si1_vs_si2_points(traces[traces.mode .== "postsynaptic", :]),
    )
    mkpath(BIO_OUTPUT_DIR)
    write_calibration_provenance()
    CSV.write(BIO_DRIVEN_TRACES_CSV, traces)
    CSV.write(BIO_DRIVEN_POINTS_CSV, points)
    println("Saved $(BIO_DRIVEN_TRACES_CSV)")
    println("Saved $(BIO_DRIVEN_POINTS_CSV)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
