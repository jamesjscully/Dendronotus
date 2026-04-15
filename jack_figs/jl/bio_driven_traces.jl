using CSV
using DataFrames
using OrdinaryDiffEq

include(joinpath(@__DIR__, "legacy", "compare_synaptic_plasticity_models.jl"))

const BIO_OUTPUT_DIR = joinpath(@__DIR__, "outputs")
const AKIRA_DIR = joinpath(dirname(@__DIR__), "akira", "extracted")
const BIO_SI1_TXT = joinpath(AKIRA_DIR, "130618-03 Si1.txt")
const BIO_SI2_TXT = joinpath(AKIRA_DIR, "130618-03 Si2.txt")
const BIO_DRIVEN_TRACES_CSV = joinpath(BIO_OUTPUT_DIR, "bio_driven_traces.csv")
const BIO_DRIVEN_POINTS_CSV = joinpath(BIO_OUTPUT_DIR, "bio_driven_points.csv")

const TRACE_PARAMS = updated_params(
    default_params();
    direct_post_base_g = 0.00064,
)

const PRESYNAPTIC_TRACE_GAIN = 3.5869003160011492
const POSTSYNAPTIC_TRACE_GAIN = 1.0303030303030303
const TRACE_SAVEAT_MS = 10.0

const REPRESENTATIVE_TRACE_U0_PRESYNAPTIC = Float64[
    -44.0, -43.29865028, -43.29889795, -43.3004827, -43.30074954,
    0.9, 0.83274431, 0.83273889, 0.83270632, 0.83270047,
    0.3, 1.12048858, 1.12049914, 1.12043197, 1.12044273,
    0.0, 0.54849525, 0.54850651, 0.54857397, 0.54858613,
    0.0, 0.09880536, 0.09880305, 0.09878894, 0.09878646,
    1.7e-7, 0.0,
    0.0, 0.00010008, 0.0001, 0.0001, 0.0001,
    0.0, 0.0, 0.0, 0.0,
    0.06775411,
]

const REPRESENTATIVE_TRACE_U0_POSTSYNAPTIC = Float64[
    -44.0, -69.70334268, -41.27715337, -48.60074388, -70.74760465,
    0.9, 0.08516567, 0.84524093, 0.80227324, 0.07747315,
    0.3, 0.78928039, 0.66296361, 0.68454718, 0.72991149,
    0.0, 0.99310944, 0.37556175, 0.32905197, 0.99397798,
    0.0, 0.00585579, 0.22486165, 0.32778749, 0.00540264,
    1.0, 0.0,
    0.0, 0.00561382, 0.00149859, 0.89754837, 0.13323279,
    3.493e-5, 0.23438674, 0.34682234, 0.00020042,
    0.06775411,
]

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
    du[2] = 4.0 * m1^3 * h1 * (30.0 - V1) + 0.3 * n1^4 * (-75.0 - V1) + 0.01 * x1 * (30.0 - V1) + 0.03 * Ca1 / (0.5 + Ca1) * (-75.0 - V1) + 0.003 * (-40.0 - V1) - g41_eff * (V1 - 30.0) * s4 - params.g12 * (V1 + 80.0) * s21 / params.scale2 + params.gelec * (V4 - V1)
    du[3] = 4.0 * m2^3 * h2 * (30.0 - V2) + 0.3 * n2^4 * (-75.0 - V2) + 0.01 * x2 * (30.0 - V2) + 0.03 * Ca2 / (0.5 + Ca2) * (-75.0 - V2) + 0.003 * (-40.0 - V2) - g32_eff * (V2 - 30.0) * s3 - params.g12 * (V2 + 80.0) * s12 / params.scale2 + params.gelec * (V3 - V2)
    du[4] = 4.0 * m3^3 * h3 * (30.0 - V3) + 0.3 * n3^4 * (-75.0 - V3) + 0.01 * x3 * (30.0 - V3) + 0.03 * Ca3 / (0.5 + Ca3) * (-75.0 - V3) + 0.003 * (-40.0 - V3) - params.g23 * (V3 + 80.0) * s2 / params.scalei - params.g34 * (V3 + 80.0) * s43 / params.scale3 + params.gelec * (V2 - V3) - params.g0 * (V3 - 30.0) * s0 / params.scale1
    du[5] = 4.0 * m4^3 * h4 * (30.0 - V4) + 0.3 * n4^4 * (-75.0 - V4) + 0.01 * x4 * (30.0 - V4) + 0.03 * Ca4 / (0.5 + Ca4) * (-75.0 - V4) + 0.003 * (-40.0 - V4) - params.g14 * (V4 + 80.0) * s1 / params.scalei - params.g34 * (V4 + 80.0) * s34 / params.scale3 + params.gelec * (V1 - V4) - params.g0 * (V4 - 30.0) * s0 / params.scale1

    du[6] = 0.0
    du[7] = ((1.0 / (exp(0.15 * (-V1 - 50.0 + params.x_shift)) + 1.0)) - x1) / 100.0
    du[8] = ((1.0 / (exp(0.15 * (-V2 - 50.0 + params.x_shift)) + 1.0)) - x2) / 100.0
    du[9] = ((1.0 / (exp(0.15 * (-V3 - 50.0 + params.x_shift)) + 1.0)) - x3) / 100.0
    du[10] = ((1.0 / (exp(0.15 * (-V4 - 50.0 + params.x_shift)) + 1.0)) - x4) / 100.0

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

    du[28] = params.alpha1 * (1.0 - s0) * syn_activation(V0) - params.beta1 * s0
    du[29] = params.alphai * s1 * (1.0 - s1) * syn_activation(V1) - params.betai * (s1 - 0.0001)
    du[30] = params.alphai * s2 * (1.0 - s2) * syn_activation(V2) - params.betai * (s2 - 0.0001)
    du[31] = alpha_exc_eff * s3 * (1.0 - s3) * syn_activation(V3) - params.betax * (s3 - 0.0001)
    du[32] = alpha_exc_eff * s4 * (1.0 - s4) * syn_activation(V4) - params.betax * (s4 - 0.0001)
    du[33] = params.alpha2 * (1.0 - s12) * syn_activation(V1) - params.beta2 * s12
    du[34] = params.alpha2 * (1.0 - s21) * syn_activation(V2) - params.beta2 * s21
    du[35] = params.alpha3 * (1.0 - s34) * syn_activation(V3) - params.beta3 * s34
    du[36] = params.alpha3 * (1.0 - s43) * syn_activation(V4) - params.beta3 * s43
    du[37] = params.alpham * sm * (1.0 - sm) * syn_activation(V0) - params.betam * (sm - 0.0001)
    return nothing
end

function representative_trace_initial_state(mode::ControlMode)
    mode == presynaptic && return copy(REPRESENTATIVE_TRACE_U0_PRESYNAPTIC)
    mode == postsynaptic && return copy(REPRESENTATIVE_TRACE_U0_POSTSYNAPTIC)
    error("Unsupported mode $(mode)")
end

function simulate_bio_driven_trace(params::ModelParams, gain::Float64, mode::ControlMode, drive_time_s::Vector{Float64}, drive_voltage_mv::Vector{Float64})
    sim_params = mode_params(params, mode)
    problem = ODEProblem(
        publication_driven_network_ode!,
        representative_trace_initial_state(mode),
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
    for i in 2:length(time_s)
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

function detect_bursts_from_voltage(time_s::AbstractVector{<:Real}, voltage_mv::AbstractVector{<:Real}; threshold_mv::Float64 = -40.0, min_down_time_s::Float64 = 1.0)
    above = Float64.(voltage_mv) .> threshold_mv
    up_idx = findall(diff(Int.(above)) .== 1) .+ 1
    down_idx = findall(diff(Int.(above)) .== -1) .+ 1
    starts = Int[]
    ends = Int[]
    for u in up_idx
        prev_down = filter(d -> d < u, down_idx)
        if !isempty(prev_down)
            d = prev_down[end]
            if Float64(time_s[u]) - Float64(time_s[d]) >= min_down_time_s
                push!(starts, u)
                next_down = filter(dn -> dn > u, down_idx)
                push!(ends, isempty(next_down) ? length(time_s) : next_down[1])
            end
        elseif Float64(time_s[u]) >= min_down_time_s
            push!(starts, u)
            next_down = filter(dn -> dn > u, down_idx)
            push!(ends, isempty(next_down) ? length(time_s) : next_down[1])
        end
    end
    return starts, ends
end

function matched_si1_vs_si2_points(trace::DataFrame)
    time_s = trace[!, :time_s]
    v0 = trace[!, :V0]
    v1 = trace[!, :V1]
    burst_starts, _ = detect_bursts_from_voltage(time_s, v1; threshold_mv = -40.0, min_down_time_s = 1.0)
    burst_times = time_s[burst_starts]
    if length(burst_times) < 2
        return DataFrame(mode = String[], si1_frequency_hz = Float64[], burst_frequency_hz = Float64[])
    end

    spike_times = threshold_crossing_spikes(time_s, v0; threshold_mv = -20.0, refractory_s = 0.05)
    if length(spike_times) < 2
        return DataFrame(mode = String[], si1_frequency_hz = Float64[], burst_frequency_hz = Float64[])
    end

    spike_freq = 1.0 ./ diff(spike_times)
    spike_freq_times = spike_times[2:end]
    burst_freq = 1.0 ./ diff(burst_times)

    si1_at_bursts = Float64[]
    for bt in burst_times[2:end]
        idx = argmin(abs.(spike_freq_times .- bt))
        push!(si1_at_bursts, spike_freq[idx])
    end

    return DataFrame(
        mode = fill(String(trace[1, :mode]), length(burst_freq)),
        si1_frequency_hz = si1_at_bursts,
        burst_frequency_hz = burst_freq,
    )
end

function main()
    bio_t1, bio_v1 = load_two_column_txt(BIO_SI1_TXT)
    bio_t1 .-= first(bio_t1)
    traces = vcat(
        simulate_bio_driven_trace(TRACE_PARAMS, PRESYNAPTIC_TRACE_GAIN, presynaptic, bio_t1, bio_v1),
        simulate_bio_driven_trace(TRACE_PARAMS, POSTSYNAPTIC_TRACE_GAIN, postsynaptic, bio_t1, bio_v1),
    )
    points = vcat(
        matched_si1_vs_si2_points(traces[traces.mode .== "presynaptic", :]),
        matched_si1_vs_si2_points(traces[traces.mode .== "postsynaptic", :]),
    )
    mkpath(BIO_OUTPUT_DIR)
    CSV.write(BIO_DRIVEN_TRACES_CSV, traces)
    CSV.write(BIO_DRIVEN_POINTS_CSV, points)
    println("Saved $(BIO_DRIVEN_TRACES_CSV)")
    println("Saved $(BIO_DRIVEN_POINTS_CSV)")
end

main()
