using DataFrames

@enum NeuromodMechanism presynaptic_kinetics postsynaptic_conductance

Base.@kwdef struct OldMatlabNeuromodParams
    mechanism::NeuromodMechanism = presynaptic_kinetics
    tmax_s::Float64 = 200.0
    step_ms::Float64 = 0.1
    sample_ms::Float64 = 2.0
    gelec::Float64 = 0.002
    alpha1::Float64 = 0.01
    beta1::Float64 = 0.002
    g0::Float64 = 0.001
    alpham::Float64 = 0.005
    betam::Float64 = 0.0001
    gm::Float64 = 3.0
    alphax::Float64 = 0.011
    betax::Float64 = 0.001
    g41::Float64 = 0.001
    g32::Float64 = 0.001
    g14::Float64 = 0.005
    g23::Float64 = 0.005
    alphai::Float64 = 0.015
    betai::Float64 = 0.0005
    alpha3::Float64 = 0.01
    beta3::Float64 = 0.002
    g34::Float64 = 0.005
    g43::Float64 = 0.005
    alpha2::Float64 = 0.01
    beta2::Float64 = 0.01
    g21::Float64 = 0.01
    g12::Float64 = 0.01
    ca_shift0::Float64 = -10.0
    ca_shift1::Float64 = -25.0
    ca_shift2::Float64 = -25.0
    ca_shift3::Float64 = -25.0
    ca_shift4::Float64 = -25.0
    x_shift::Float64 = -4.0
    iapp::Float64 = 0.0
    t1_ms::Float64 = 1000.0
    t2_ms::Float64 = 35000.0
    gh::Float64 = 0.0005
    vhh::Float64 = -53.0
    v0_init::Float64 = -44.0
    v1_init::Float64 = -44.0
    v2_init::Float64 = -44.0
    v3_init::Float64 = -44.0
    v4_init::Float64 = -44.0
    x0_init::Float64 = 0.9
    x1_init::Float64 = 0.6
    x2_init::Float64 = 0.6
    x3_init::Float64 = 0.5
    x4_init::Float64 = 0.6
    ca0_init::Float64 = 0.3
    ca1_init::Float64 = 1.0
    ca2_init::Float64 = 1.0
    ca3_init::Float64 = 1.1
    ca4_init::Float64 = 1.0
    s34_init::Float64 = 0.1
    s43_init::Float64 = 0.1
    sm_init::Float64 = 0.8
end

function old_matlab_presynaptic_params(; kwargs...)
    OldMatlabNeuromodParams(;
        mechanism = presynaptic_kinetics,
        g0 = 0.001,
        g41 = 0.001,
        g32 = 0.001,
        x0_init = 0.9,
        sm_init = 0.8,
        kwargs...,
    )
end

function old_matlab_postsynaptic_params(; kwargs...)
    OldMatlabNeuromodParams(;
        mechanism = postsynaptic_conductance,
        g0 = 0.002,
        g41 = 0.04,
        g32 = 0.04,
        x0_init = 0.0,
        sm_init = 0.99,
        kwargs...,
    )
end

scale1(p::OldMatlabNeuromodParams) = p.alpha1 / (p.alpha1 + p.beta1)
scalem(p::OldMatlabNeuromodParams) = (p.alpham - p.betam) / p.alpham
scalex(p::OldMatlabNeuromodParams) = (p.alphax - p.betax) / p.alphax
scalei(p::OldMatlabNeuromodParams) = (p.alphai - p.betai) / p.alphai
scale2(p::OldMatlabNeuromodParams) = p.alpha2 / (p.alpha2 + p.beta2)
scale3(p::OldMatlabNeuromodParams) = p.alpha3 / (p.alpha3 + p.beta3)

heaviside(x) = ifelse(x >= 0, 1.0, 0.0)
sigma_syn(v) = 1.0 / (1.0 + exp(-20.0 * (v + 20.0)))

function alpha_m_gate(v)
    u = 127.0 * v / 105.0 + 8265.0 / 105.0
    z = 50.0 - u
    abs(z) < 1e-8 && return 1.0
    return 0.1 * z / (exp(z / 10.0) - 1.0)
end

function sodium_activation(v)
    am = alpha_m_gate(v)
    bm = 4.0 * exp((25.0 - (127.0 * v / 105.0 + 8265.0 / 105.0)) / 18.0)
    return am / (am + bm)
end

function alpha_n_gate(v)
    u = 127.0 * v / 105.0 + 8265.0 / 105.0
    z = 55.0 - u
    abs(z) < 1e-8 && return 0.1
    return 0.01 * z / (exp(z / 10.0) - 1.0)
end

function intrinsic_voltage_rhs(v, h, n, x, ca)
    m = sodium_activation(v)
    4.0 * m^3 * h * (30.0 - v) +
        0.3 * n^4 * (-75.0 - v) +
        0.01 * x * (30.0 - v) +
        0.03 * ca / (0.5 + ca) * (-75.0 - v) +
        0.003 * (-40.0 - v)
end

function calcium_rhs(v, x, ca, ca_shift, rate)
    rate * (0.0085 * x * (140.0 - v + ca_shift) - ca)
end

x_rhs(v, x, x_shift) = ((1.0 / (exp(0.15 * (-v - 50.0 + x_shift)) + 1.0)) - x) / 100.0

function h_rhs(v, h)
    u = 127.0 * v / 105.0 + 8265.0 / 105.0
    ((1.0 - h) * (0.07 * exp((25.0 - u) / 20.0)) - h * (1.0 / (1.0 + exp((55.0 - u) / 10.0)))) / 12.5
end

function n_rhs(v, n)
    u = 127.0 * v / 105.0 + 8265.0 / 105.0
    ((1.0 - n) * alpha_n_gate(v) - n * (0.125 * exp((45.0 - u) / 80.0))) / 12.5
end

y_rhs(v, y, vhh) = 0.5 * ((1.0 / (1.0 + exp(10.0 * (v - vhh)))) - y) / (7.1 + 10.4 / (1.0 + exp((v + 68.0) / 2.2)))

function simulate_old_matlab_neuromod(p::OldMatlabNeuromodParams)
    nt = round(Int, p.tmax_s * 1000.0 / p.step_ms)
    sample_stride = round(Int, p.sample_ms / p.step_ms)
    nt_save = div(nt, sample_stride)

    time_ms = Vector{Float64}(undef, nt_save)
    V0s = similar(time_ms); V1s = similar(time_ms); V2s = similar(time_ms); V3s = similar(time_ms); V4s = similar(time_ms)
    s0s = similar(time_ms); sms = similar(time_ms); s1s = similar(time_ms); s2s = similar(time_ms); s3s = similar(time_ms); s4s = similar(time_ms)
    s12s = similar(time_ms); s21s = similar(time_ms); s34s = similar(time_ms); s43s = similar(time_ms)
    ca1s = similar(time_ms); ca2s = similar(time_ms); ca3s = similar(time_ms); ca4s = similar(time_ms)

    V0 = p.v0_init; V1 = p.v1_init; V2 = p.v2_init; V3 = p.v3_init; V4 = p.v4_init
    x0 = p.x0_init; x1 = p.x1_init; x2 = p.x2_init; x3 = p.x3_init; x4 = p.x4_init
    Ca0 = p.ca0_init; Ca1 = p.ca1_init; Ca2 = p.ca2_init; Ca3 = p.ca3_init; Ca4 = p.ca4_init
    h0 = 0.0; h1 = 0.0; h2 = 0.0; h3 = 0.0; h4 = 0.0
    n0 = 0.0; n1 = 0.0; n2 = 0.0; n3 = 0.0; n4 = 0.0
    y1 = 0.0; y2 = 0.0
    s0 = 0.0; sm = p.sm_init; s1 = 0.0; s2 = 0.0; s3 = 0.0; s4 = 0.0
    s12 = 0.0; s21 = 0.0; s34 = p.s34_init; s43 = p.s43_init

    sc1 = scale1(p); scm = scalem(p); scx = scalex(p); sci = scalei(p); sc2 = scale2(p); sc3 = scale3(p)

    j = 0
    for i in 1:nt
        tt = (i - 1) * p.step_ms

        drive = p.iapp * heaviside(tt - p.t1_ms) * heaviside(p.t2_ms + 1250.0 - tt)
        V0new = V0 + p.step_ms * (drive + intrinsic_voltage_rhs(V0, h0, n0, x0, Ca0))

        if p.mechanism == presynaptic_kinetics
            exc_factor = 1.0
            exc_scale = 1.0
        else
            exc_factor = 1.0 + p.gm * sm / scm
            exc_scale = scx
        end

        V1new = V1 + p.step_ms * (intrinsic_voltage_rhs(V1, h1, n1, x1, Ca1) -
            p.g41 * (V1 - 30.0) * s4 * exc_factor / exc_scale -
            p.g12 * (V1 + 80.0) * s21 / sc2 +
            p.gelec * (V4 - V1))
        V2new = V2 + p.step_ms * (intrinsic_voltage_rhs(V2, h2, n2, x2, Ca2) -
            p.g32 * (V2 - 30.0) * s3 * exc_factor / exc_scale -
            p.g12 * (V2 + 80.0) * s12 / sc2 +
            p.gelec * (V3 - V2))
        V3new = V3 + p.step_ms * (intrinsic_voltage_rhs(V3, h3, n3, x3, Ca3) -
            p.g23 * (V3 + 80.0) * s2 / sci -
            p.g34 * (V3 + 80.0) * s43 / sc3 +
            p.gelec * (V2 - V3) -
            p.g0 * (V3 - 30.0) * s0 / sc1)
        V4new = V4 + p.step_ms * (intrinsic_voltage_rhs(V4, h4, n4, x4, Ca4) -
            p.g14 * (V4 + 80.0) * s1 / sci -
            p.g34 * (V4 + 80.0) * s34 / sc3 +
            p.gelec * (V1 - V4) -
            p.g0 * (V4 - 30.0) * s0 / sc1)

        Ca0new = Ca0 + p.step_ms * calcium_rhs(V0, x0, Ca0, p.ca_shift0, 0.000012)
        Ca1new = Ca1 + p.step_ms * calcium_rhs(V1, x1, Ca1, p.ca_shift1, 0.0003)
        Ca2new = Ca2 + p.step_ms * calcium_rhs(V2, x2, Ca2, p.ca_shift2, 0.0003)
        Ca3new = Ca3 + p.step_ms * calcium_rhs(V3, x3, Ca3, p.ca_shift3, 0.0003)
        Ca4new = Ca4 + p.step_ms * calcium_rhs(V4, x4, Ca4, p.ca_shift4, 0.0003)

        x0new = x0 + p.step_ms * x_rhs(V0, x0, p.x_shift)
        x1new = x1 + p.step_ms * x_rhs(V1, x1, p.x_shift)
        x2new = x2 + p.step_ms * x_rhs(V2, x2, p.x_shift)
        x3new = x3 + p.step_ms * x_rhs(V3, x3, p.x_shift)
        x4new = x4 + p.step_ms * x_rhs(V4, x4, p.x_shift)

        h0new = h0 + p.step_ms * h_rhs(V0, h0)
        h1new = h1 + p.step_ms * h_rhs(V1, h1)
        h2new = h2 + p.step_ms * h_rhs(V2, h2)
        h3new = h3 + p.step_ms * h_rhs(V3, h3)
        h4new = h4 + p.step_ms * h_rhs(V4, h4)

        n0new = n0 + p.step_ms * n_rhs(V0, n0)
        n1new = n1 + p.step_ms * n_rhs(V1, n1)
        n2new = n2 + p.step_ms * n_rhs(V2, n2)
        n3new = n3 + p.step_ms * n_rhs(V3, n3)
        n4new = n4 + p.step_ms * n_rhs(V4, n4)

        s0new = s0 + p.step_ms * (p.alpha1 * (1.0 - s0) * sigma_syn(V0) - p.beta1 * s0)
        smnew = sm + p.step_ms * (p.alpham * sm * (1.0 - sm) * sigma_syn(V0) - p.betam * (sm - 0.0001))

        if p.mechanism == presynaptic_kinetics
            alphax_eff = p.alphax * (1.0 + p.gm * sm) / scm
        else
            alphax_eff = p.alphax
        end
        s4new = s4 + p.step_ms * (alphax_eff * s4 * (1.0 - s4) * sigma_syn(V4) - p.betax * (s4 - 0.0001))
        s3new = s3 + p.step_ms * (alphax_eff * s3 * (1.0 - s3) * sigma_syn(V3) - p.betax * (s3 - 0.0001))

        s1new = s1 + p.step_ms * (p.alphai * s1 * (1.0 - s1) * sigma_syn(V1) - p.betai * (s1 - 0.0001))
        s2new = s2 + p.step_ms * (p.alphai * s2 * (1.0 - s2) * sigma_syn(V2) - p.betai * (s2 - 0.0001))
        s12new = s12 + p.step_ms * (p.alpha2 * (1.0 - s12) * sigma_syn(V1) - p.beta2 * s12)
        s21new = s21 + p.step_ms * (p.alpha2 * (1.0 - s21) * sigma_syn(V2) - p.beta2 * s21)
        s34new = s34 + p.step_ms * (p.alpha3 * (1.0 - s34) * sigma_syn(V3) - p.beta3 * s34)
        s43new = s43 + p.step_ms * (p.alpha3 * (1.0 - s43) * sigma_syn(V4) - p.beta3 * s43)

        y1new = y1 + p.step_ms * y_rhs(V1, y1, p.vhh)
        y2new = y2 + p.step_ms * y_rhs(V2, y2, p.vhh)

        V0 = V0new; V1 = V1new; V2 = V2new; V3 = V3new; V4 = V4new
        x0 = x0new; x1 = x1new; x2 = x2new; x3 = x3new; x4 = x4new
        Ca0 = Ca0new; Ca1 = Ca1new; Ca2 = Ca2new; Ca3 = Ca3new; Ca4 = Ca4new
        h0 = h0new; h1 = h1new; h2 = h2new; h3 = h3new; h4 = h4new
        n0 = n0new; n1 = n1new; n2 = n2new; n3 = n3new; n4 = n4new
        y1 = y1new; y2 = y2new
        s0 = s0new; sm = smnew; s1 = s1new; s2 = s2new; s3 = s3new; s4 = s4new
        s12 = s12new; s21 = s21new; s34 = s34new; s43 = s43new

        if i % sample_stride == 0
            j += 1
            time_ms[j] = tt
            V0s[j] = V0; V1s[j] = V1; V2s[j] = V2; V3s[j] = V3; V4s[j] = V4
            s0s[j] = s0; sms[j] = sm; s1s[j] = s1; s2s[j] = s2; s3s[j] = s3; s4s[j] = s4
            s12s[j] = s12; s21s[j] = s21; s34s[j] = s34; s43s[j] = s43
            ca1s[j] = Ca1; ca2s[j] = Ca2; ca3s[j] = Ca3; ca4s[j] = Ca4
        end
    end

    DataFrame(
        time_s = time_ms ./ 1000.0,
        V0 = V0s, V1 = V1s, V2 = V2s, V3 = V3s, V4 = V4s,
        s0 = s0s, sm = sms, s1 = s1s, s2 = s2s, s3 = s3s, s4 = s4s,
        s12 = s12s, s21 = s21s, s34 = s34s, s43 = s43s,
        Ca1 = ca1s, Ca2 = ca2s, Ca3 = ca3s, Ca4 = ca4s,
    )
end

function threshold_crossings(time_s, voltage_mv; threshold_mv = -20.0, refractory_s = 0.02)
    spikes = Float64[]
    last_spike = -Inf
    above_prev = voltage_mv[1] >= threshold_mv
    for i in 2:length(voltage_mv)
        above = voltage_mv[i] >= threshold_mv
        if above && !above_prev && time_s[i] - last_spike >= refractory_s
            push!(spikes, time_s[i])
            last_spike = time_s[i]
        end
        above_prev = above
    end
    return spikes
end

function burst_start_times_from_spikes(spikes::Vector{Float64}; gap_s = 0.45, min_spikes = 3)
    isempty(spikes) && return Float64[]
    starts = Float64[]
    run_start = spikes[1]
    run_count = 1
    for i in 2:length(spikes)
        if spikes[i] - spikes[i - 1] <= gap_s
            run_count += 1
        else
            run_count >= min_spikes && push!(starts, run_start)
            run_start = spikes[i]
            run_count = 1
        end
    end
    run_count >= min_spikes && push!(starts, run_start)
    return starts
end

function trace_rate_summary(trace::DataFrame; si2_column = :V1, si1_column = :V0, t_start_s = 50.0, t_stop_s = 190.0)
    time = trace.time_s
    keep = findall(t -> t_start_s <= t <= t_stop_s, time)
    isempty(keep) && return (si1_frequency_hz = NaN, burst_frequency_hz = NaN, n_bursts = 0)
    si1_spikes = threshold_crossings(time[keep], trace[keep, si1_column])
    bursts = threshold_crossings(time[keep], trace[keep, si2_column]; refractory_s = 2.5)
    burst_frequency = length(bursts) >= 2 ? (length(bursts) - 1) / (last(bursts) - first(bursts)) : NaN
    si1_frequency = length(si1_spikes) / (time[keep[end]] - time[keep[1]])
    return (si1_frequency_hz = si1_frequency, burst_frequency_hz = burst_frequency, n_bursts = length(bursts))
end
