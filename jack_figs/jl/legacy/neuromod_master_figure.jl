using CSV
using CairoMakie
using DataFrames
using FileIO
using GLMakie
using Statistics
using OrdinaryDiffEq

include(joinpath(@__DIR__, "publication_plasticity_scan.jl"))

const MASTER_PNG = joinpath(OUTPUT_DIR, "neuromod_master_figure.png")
const MASTER_MODEL_POINTS_CSV = joinpath(OUTPUT_DIR, "neuromod_master_model_points.csv")
const MASTER_MODEL_TRACES_CSV = joinpath(OUTPUT_DIR, "neuromod_master_model_traces.csv")
const MASTER_DRIVEN_TRACES_CSV = joinpath(OUTPUT_DIR, "neuromod_master_driven_traces.csv")
const PAPER_FIG5_JPG = joinpath(dirname(dirname(@__DIR__)), "paper", "fig5.jpg")
const LEGACY_CIRCUIT_SOURCE = joinpath(dirname(@__DIR__), "s1_Fig2A_drives_cpg_als_circuit.jpg")
const LEGACY_CIRCUIT_CROP = joinpath(OUTPUT_DIR, "legacy_circuit_crop2.jpg")

const AKIRA_DIR = joinpath(dirname(@__DIR__), "akira")
const AKIRA_EXTRACTED = joinpath(AKIRA_DIR, "extracted")
const BIO_SI1_TXT = joinpath(AKIRA_EXTRACTED, "130618-03 Si1.txt")
const BIO_SI2_TXT = joinpath(AKIRA_EXTRACTED, "130618-03 Si2.txt")
const BIO_FIG2C_TSV = joinpath(AKIRA_DIR, "xlsx_cache", "Fig_2C.tsv")
const DRIVEN_TRACE_SAVEAT_MS = 10.0
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

const MODE_COLORS = Dict(
    "presynaptic" => RGBf(0.12, 0.42, 0.86),
    "postsynaptic" => RGBf(0.80, 0.22, 0.18),
)

const MODE_MARKERS = Dict(
    "presynaptic" => :circle,
    "postsynaptic" => :rect,
)

struct CompareModelConfig
    mode_name::String
    g0::Float64
    g41::Float64
    g32::Float64
    sm_g0::Float64
    sm_a0::Float64
    x0_init::Float64
    use_g_mod::Bool
    use_a_mod::Bool
end

function compare_model_config(mode_name::String)
    if mode_name == "postsynaptic"
        return CompareModelConfig(mode_name, 0.002, 0.04, 0.04, 0.99, 0.0, 0.9, true, false)
    elseif mode_name == "presynaptic"
        return CompareModelConfig(mode_name, 0.001, 0.001, 0.001, 0.0, 0.8, 0.9, false, true)
    end
    error("Unsupported compare-model mode $(mode_name)")
end

clean_name(x) = replace(String(x), "_" => " ")

function finite_trace_cache(df::DataFrame)
    return all(isfinite, Float64.(df[!, :V0])) && all(isfinite, Float64.(df[!, :V1]))
end

function ensure_publication_cache()
    if !isfile(PUB_SCAN_SUMMARY_CSV) || !isfile(PUB_SCAN_RAW_CSV)
        raw, summary = generate_publication_data()
        save_publication_data(raw, summary)
    end
end

function publication_summary_subset()
    summary = load_publication_summary()
    return summary[in.(summary[!, :mode], Ref(["presynaptic", "postsynaptic"])), :]
end

function publication_raw_subset()
    raw = load_publication_raw()
    return raw[in.(raw[!, :mode], Ref(["presynaptic", "postsynaptic"])), :]
end

function representative_gain(summary::DataFrame, mode_name::String; target_hz::Float64 = 0.12)
    subset = summary[(summary[!, :mode] .== mode_name) .& (summary[!, :n_pre_intervals] .> 0), :]
    isempty(subset) && error("No valid publication-summary points for $(mode_name)")
    idx = argmin(abs.(Float64.(subset[!, :mean_pre_frequency_hz]) .- target_hz))
    return Float64(subset[idx, :control_gain])
end

function simulation_dataframe(sim, mode_name::String)
    return DataFrame(
        mode = fill(mode_name, length(sim.time_s)),
        time_s = sim.time_s,
        V0 = sim.V0,
        V1 = sim.V1,
    )
end

function finite_xy(x::AbstractVector, y::AbstractVector)
    keep = isfinite.(Float64.(x)) .& isfinite.(Float64.(y))
    return Float64.(x[keep]), Float64.(y[keep])
end

function load_two_column_txt(path::AbstractString)
    df = DataFrame(CSV.File(path; delim = '\t', normalizenames = false))
    x = Float64.(df[!, 1])
    y = Float64.(df[!, 2])
    return x, y
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

function upward_spike_times(time_s::AbstractVector{<:Real}, voltage_mv::AbstractVector{<:Real}; threshold_mv::Float64, refractory_s::Float64)
    spikes = Float64[]
    last_spike = -Inf
    for i in 2:length(time_s)
        v_prev = Float64(voltage_mv[i - 1])
        v_now = Float64(voltage_mv[i])
        if v_prev < threshold_mv && v_now >= threshold_mv
            t_prev = Float64(time_s[i - 1])
            t_now = Float64(time_s[i])
            frac = (threshold_mv - v_prev) / max(v_now - v_prev, eps())
            t_cross = t_prev + frac * (t_now - t_prev)
            if isempty(spikes) || (t_cross - last_spike) >= refractory_s
                push!(spikes, t_cross)
                last_spike = t_cross
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

function representative_burst_window(time_s::AbstractVector{<:Real}, voltage_mv::AbstractVector{<:Real}; threshold_mv::Float64 = -40.0, min_down_time_s::Float64 = 1.0, pad_s::Float64 = 1.0)
    starts, ends = detect_bursts_from_voltage(time_s, voltage_mv; threshold_mv = threshold_mv, min_down_time_s = min_down_time_s)
    isempty(starts) && return (Float64(first(time_s)), Float64(last(time_s)))
    idx = clamp(cld(length(starts), 2), 1, length(starts))
    t_start = max(Float64(first(time_s)), Float64(time_s[starts[idx]]) - pad_s)
    t_stop = min(Float64(last(time_s)), Float64(time_s[ends[idx]]) + pad_s)
    return t_start, t_stop
end

function first_burst_onset(time_s::AbstractVector{<:Real}, voltage_mv::AbstractVector{<:Real}; threshold_mv::Float64 = -40.0, min_down_time_s::Float64 = 1.0)
    starts, _ = detect_bursts_from_voltage(time_s, voltage_mv; threshold_mv = threshold_mv, min_down_time_s = min_down_time_s)
    isempty(starts) && return NaN
    return Float64(time_s[first(starts)])
end

function simulate_compare_model(mode_name::String; tmax_s::Float64 = 200.0, step_ms::Float64 = 0.1, saveat_ms::Float64 = 2.0)
    cfg = compare_model_config(mode_name)

    alpha1 = 0.01
    beta1 = 0.002
    scale1 = alpha1 / (alpha1 + beta1)
    alpham_g = 0.005
    betam_g = 0.0001
    scalem_g = (alpham_g - betam_g) / alpham_g
    gm = 3.0
    alpham_a = 0.005
    betam_a = 0.0001
    scalem_a = (alpham_a - betam_a) / alpham_a
    gm_a = 3.0
    alphax = 0.011
    betax = 0.001
    alphai = 0.015
    betai = 0.0005
    scalei = (alphai - betai) / alphai
    alpha3 = 0.01
    beta3 = 0.002
    scale3 = alpha3 / (alpha3 + beta3)
    alpha2 = 0.01
    beta2 = 0.01
    scale2 = alpha2 / (alpha2 + beta2)
    gelec = 0.002
    g14 = 0.005
    g23 = 0.005
    g34 = 0.005
    g43 = 0.005
    g21 = 0.01
    g12 = 0.01
    Ca_shift0 = -10.0
    Ca_shift4 = -25.0
    Ca_shift1 = -25.0
    Ca_shift3 = Ca_shift4
    Ca_shift2 = Ca_shift1
    x_shift = -4.0
    Vhh = -53.0

    V0 = -44.0
    V1 = -44.0
    V2 = -44.0
    V3 = -44.0
    V4 = -44.0
    x0 = cfg.x0_init
    x1 = 0.6
    x2 = 0.6
    x3 = 0.5
    x4 = 0.6
    Ca0 = 0.3
    Ca1 = 1.0
    Ca2 = 1.0
    Ca3 = 1.1
    Ca4 = 1.0
    h0 = 0.0
    h1 = 0.0
    h2 = 0.0
    h3 = 0.0
    h4 = 0.0
    n0 = 0.0
    n1 = 0.0
    n2 = 0.0
    n3 = 0.0
    n4 = 0.0
    y1 = 0.0
    y2 = 0.0
    s0 = 0.0
    s1 = 0.0
    s2 = 0.0
    s3 = 0.0
    s4 = 0.0
    s12 = 0.0
    s21 = 0.0
    s34 = 0.1
    s43 = 0.1
    sm_g = cfg.sm_g0
    sm_a = cfg.sm_a0

    nt = round(Int, tmax_s * 1000 / step_ms)
    isave = round(Int, saveat_ms / step_ms)
    nsave = round(Int, tmax_s * 1000 / saveat_ms) + 1

    time_s = zeros(Float64, nsave)
    vv0 = zeros(Float64, nsave)
    vv1 = zeros(Float64, nsave)
    vv2 = zeros(Float64, nsave)
    vv3 = zeros(Float64, nsave)
    vv4 = zeros(Float64, nsave)
    save_idx = 1
    time_s[1] = 0.0
    vv0[1] = V0
    vv1[1] = V1
    vv2[1] = V2
    vv3[1] = V3
    vv4[1] = V4

    for i in 1:nt
        Vs0 = 127.0 * V0 / 105.0 + 8265.0 / 105.0
        Vs1 = 127.0 * V1 / 105.0 + 8265.0 / 105.0
        Vs2 = 127.0 * V2 / 105.0 + 8265.0 / 105.0
        Vs3 = 127.0 * V3 / 105.0 + 8265.0 / 105.0
        Vs4 = 127.0 * V4 / 105.0 + 8265.0 / 105.0

        m0 = (0.1 * (50.0 - Vs0) / (exp((50.0 - Vs0) / 10.0) - 1.0)) / ((0.1 * (50.0 - Vs0) / (exp((50.0 - Vs0) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs0) / 18.0))
        m1 = (0.1 * (50.0 - Vs1) / (exp((50.0 - Vs1) / 10.0) - 1.0)) / ((0.1 * (50.0 - Vs1) / (exp((50.0 - Vs1) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs1) / 18.0))
        m2 = (0.1 * (50.0 - Vs2) / (exp((50.0 - Vs2) / 10.0) - 1.0)) / ((0.1 * (50.0 - Vs2) / (exp((50.0 - Vs2) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs2) / 18.0))
        m3 = (0.1 * (50.0 - Vs3) / (exp((50.0 - Vs3) / 10.0) - 1.0)) / ((0.1 * (50.0 - Vs3) / (exp((50.0 - Vs3) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs3) / 18.0))
        m4 = (0.1 * (50.0 - Vs4) / (exp((50.0 - Vs4) / 10.0) - 1.0)) / ((0.1 * (50.0 - Vs4) / (exp((50.0 - Vs4) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs4) / 18.0))

        g41_eff = cfg.use_g_mod ? cfg.g41 * (1.0 + gm * sm_g / scalem_g) : cfg.g41
        g32_eff = cfg.use_g_mod ? cfg.g32 * (1.0 + gm * sm_g / scalem_g) : cfg.g32
        alpha_exc = cfg.use_a_mod ? alphax * (1.0 + gm_a * sm_a / scalem_a) : alphax

        V0_new = V0 + step_ms * (
            4.0 * m0^3 * h0 * (30.0 - V0) +
            0.3 * n0^4 * (-75.0 - V0) +
            0.01 * x0 * (30.0 - V0) +
            0.03 * Ca0 / (0.5 + Ca0) * (-75.0 - V0) +
            0.003 * (-40.0 - V0)
        )
        V1_new = V1 + step_ms * (
            4.0 * m1^3 * h1 * (30.0 - V1) +
            0.3 * n1^4 * (-75.0 - V1) +
            0.01 * x1 * (30.0 - V1) +
            0.03 * Ca1 / (0.5 + Ca1) * (-75.0 - V1) +
            0.003 * (-40.0 - V1) -
            g41_eff * (V1 - 30.0) * s4 -
            g12 * (V1 + 80.0) * s21 / scale2 +
            gelec * (V4 - V1)
        )
        V2_new = V2 + step_ms * (
            4.0 * m2^3 * h2 * (30.0 - V2) +
            0.3 * n2^4 * (-75.0 - V2) +
            0.01 * x2 * (30.0 - V2) +
            0.03 * Ca2 / (0.5 + Ca2) * (-75.0 - V2) +
            0.003 * (-40.0 - V2) -
            g32_eff * (V2 - 30.0) * s3 -
            g12 * (V2 + 80.0) * s12 / scale2 +
            gelec * (V3 - V2)
        )
        V3_new = V3 + step_ms * (
            4.0 * m3^3 * h3 * (30.0 - V3) +
            0.3 * n3^4 * (-75.0 - V3) +
            0.01 * x3 * (30.0 - V3) +
            0.03 * Ca3 / (0.5 + Ca3) * (-75.0 - V3) +
            0.003 * (-40.0 - V3) -
            g23 * (V3 + 80.0) * s2 / scalei -
            g34 * (V3 + 80.0) * s43 / scale3 +
            gelec * (V2 - V3) -
            cfg.g0 * (V3 - 30.0) * s0 / scale1
        )
        V4_new = V4 + step_ms * (
            4.0 * m4^3 * h4 * (30.0 - V4) +
            0.3 * n4^4 * (-75.0 - V4) +
            0.01 * x4 * (30.0 - V4) +
            0.03 * Ca4 / (0.5 + Ca4) * (-75.0 - V4) +
            0.003 * (-40.0 - V4) -
            g14 * (V4 + 80.0) * s1 / scalei -
            g34 * (V4 + 80.0) * s34 / scale3 +
            gelec * (V1 - V4) -
            cfg.g0 * (V4 - 30.0) * s0 / scale1
        )

        Ca0_new = Ca0 + step_ms * (0.000012 * (0.0085 * x0 * (140.0 - V0 + Ca_shift0) - Ca0))
        Ca1_new = Ca1 + step_ms * (0.0003 * (0.0085 * x1 * (140.0 - V1 + Ca_shift1) - Ca1))
        Ca2_new = Ca2 + step_ms * (0.0003 * (0.0085 * x2 * (140.0 - V2 + Ca_shift2) - Ca2))
        Ca3_new = Ca3 + step_ms * (0.0003 * (0.0085 * x3 * (140.0 - V3 + Ca_shift3) - Ca3))
        Ca4_new = Ca4 + step_ms * (0.0003 * (0.0085 * x4 * (140.0 - V4 + Ca_shift4) - Ca4))

        x0_new = x0 + step_ms * (((1.0 / (exp(0.15 * (-V0 - 50.0 + x_shift)) + 1.0)) - x0) / 100.0)
        x1_new = x1 + step_ms * (((1.0 / (exp(0.15 * (-V1 - 50.0 + x_shift)) + 1.0)) - x1) / 100.0)
        x2_new = x2 + step_ms * (((1.0 / (exp(0.15 * (-V2 - 50.0 + x_shift)) + 1.0)) - x2) / 100.0)
        x3_new = x3 + step_ms * (((1.0 / (exp(0.15 * (-V3 - 50.0 + x_shift)) + 1.0)) - x3) / 100.0)
        x4_new = x4 + step_ms * (((1.0 / (exp(0.15 * (-V4 - 50.0 + x_shift)) + 1.0)) - x4) / 100.0)

        h0_new = h0 + step_ms * (((1.0 - h0) * (0.07 * exp((25.0 - Vs0) / 20.0)) - h0 * (1.0 / (1.0 + exp((55.0 - Vs0) / 10.0)))) / 12.5)
        h1_new = h1 + step_ms * (((1.0 - h1) * (0.07 * exp((25.0 - Vs1) / 20.0)) - h1 * (1.0 / (1.0 + exp((55.0 - Vs1) / 10.0)))) / 12.5)
        h2_new = h2 + step_ms * (((1.0 - h2) * (0.07 * exp((25.0 - Vs2) / 20.0)) - h2 * (1.0 / (1.0 + exp((55.0 - Vs2) / 10.0)))) / 12.5)
        h3_new = h3 + step_ms * (((1.0 - h3) * (0.07 * exp((25.0 - Vs3) / 20.0)) - h3 * (1.0 / (1.0 + exp((55.0 - Vs3) / 10.0)))) / 12.5)
        h4_new = h4 + step_ms * (((1.0 - h4) * (0.07 * exp((25.0 - Vs4) / 20.0)) - h4 * (1.0 / (1.0 + exp((55.0 - Vs4) / 10.0)))) / 12.5)

        n0_new = n0 + step_ms * (((1.0 - n0) * (0.01 * (55.0 - Vs0) / (exp((55.0 - Vs0) / 10.0) - 1.0)) - n0 * (0.125 * exp((45.0 - Vs0) / 80.0))) / 12.5)
        n1_new = n1 + step_ms * (((1.0 - n1) * (0.01 * (55.0 - Vs1) / (exp((55.0 - Vs1) / 10.0) - 1.0)) - n1 * (0.125 * exp((45.0 - Vs1) / 80.0))) / 12.5)
        n2_new = n2 + step_ms * (((1.0 - n2) * (0.01 * (55.0 - Vs2) / (exp((55.0 - Vs2) / 10.0) - 1.0)) - n2 * (0.125 * exp((45.0 - Vs2) / 80.0))) / 12.5)
        n3_new = n3 + step_ms * (((1.0 - n3) * (0.01 * (55.0 - Vs3) / (exp((55.0 - Vs3) / 10.0) - 1.0)) - n3 * (0.125 * exp((45.0 - Vs3) / 80.0))) / 12.5)
        n4_new = n4 + step_ms * (((1.0 - n4) * (0.01 * (55.0 - Vs4) / (exp((55.0 - Vs4) / 10.0) - 1.0)) - n4 * (0.125 * exp((45.0 - Vs4) / 80.0))) / 12.5)

        y1_new = y1 + step_ms * (0.5 * ((1.0 / (1.0 + exp(10.0 * (V1 - Vhh)))) - y1) / (7.1 + 10.4 / (1.0 + exp((V1 + 68.0) / 2.2))))
        y2_new = y2 + step_ms * (0.5 * ((1.0 / (1.0 + exp(10.0 * (V2 - Vhh)))) - y2) / (7.1 + 10.4 / (1.0 + exp((V2 + 68.0) / 2.2))))

        act0 = 1.0 / (1.0 + exp(-20.0 * (V0 + 20.0)))
        act1 = 1.0 / (1.0 + exp(-20.0 * (V1 + 20.0)))
        act2 = 1.0 / (1.0 + exp(-20.0 * (V2 + 20.0)))
        act3 = 1.0 / (1.0 + exp(-20.0 * (V3 + 20.0)))
        act4 = 1.0 / (1.0 + exp(-20.0 * (V4 + 20.0)))

        s0_new = s0 + step_ms * (alpha1 * (1.0 - s0) * act0 - beta1 * s0)
        sm_g_new = cfg.use_g_mod ? sm_g + step_ms * (alpham_g * sm_g * (1.0 - sm_g) * act0 - betam_g * (sm_g - 0.0001)) : sm_g
        sm_a_new = cfg.use_a_mod ? sm_a + step_ms * (alpham_a * sm_a * (1.0 - sm_a) * act0 - betam_a * (sm_a - 0.0001)) : sm_a
        s1_new = s1 + step_ms * (alphai * s1 * (1.0 - s1) * act1 - betai * (s1 - 0.0001))
        s2_new = s2 + step_ms * (alphai * s2 * (1.0 - s2) * act2 - betai * (s2 - 0.0001))
        s3_new = s3 + step_ms * (alpha_exc * s3 * (1.0 - s3) * act3 - betax * (s3 - 0.0001))
        s4_new = s4 + step_ms * (alpha_exc * s4 * (1.0 - s4) * act4 - betax * (s4 - 0.0001))
        s12_new = s12 + step_ms * (alpha2 * (1.0 - s12) * act1 - beta2 * s12)
        s21_new = s21 + step_ms * (alpha2 * (1.0 - s21) * act2 - beta2 * s21)
        s34_new = s34 + step_ms * (alpha3 * (1.0 - s34) * act3 - beta3 * s34)
        s43_new = s43 + step_ms * (alpha3 * (1.0 - s43) * act4 - beta3 * s43)

        V0, V1, V2, V3, V4 = V0_new, V1_new, V2_new, V3_new, V4_new
        x0, x1, x2, x3, x4 = x0_new, x1_new, x2_new, x3_new, x4_new
        Ca0, Ca1, Ca2, Ca3, Ca4 = Ca0_new, Ca1_new, Ca2_new, Ca3_new, Ca4_new
        h0, h1, h2, h3, h4 = h0_new, h1_new, h2_new, h3_new, h4_new
        n0, n1, n2, n3, n4 = n0_new, n1_new, n2_new, n3_new, n4_new
        y1, y2 = y1_new, y2_new
        s0, s1, s2, s3, s4 = s0_new, s1_new, s2_new, s3_new, s4_new
        s12, s21, s34, s43 = s12_new, s21_new, s34_new, s43_new
        sm_g, sm_a = sm_g_new, sm_a_new

        if i % isave == 0
            save_idx += 1
            time_s[save_idx] = i * step_ms / 1000.0
            vv0[save_idx] = V0
            vv1[save_idx] = V1
            vv2[save_idx] = V2
            vv3[save_idx] = V3
            vv4[save_idx] = V4
        end
    end

    return DataFrame(
        mode = fill(mode_name, save_idx),
        time_s = time_s[1:save_idx],
        V0 = vv0[1:save_idx],
        V1 = vv1[1:save_idx],
        V2 = vv2[1:save_idx],
        V3 = vv3[1:save_idx],
        V4 = vv4[1:save_idx],
    )
end

function interpolate_drive(time_s::Float64, drive_t::Vector{Float64}, drive_v::Vector{Float64})
    time_s <= drive_t[1] && return drive_v[1]
    time_s >= drive_t[end] && return drive_v[end]
    idx = searchsortedlast(drive_t, time_s)
    idx = clamp(idx, 1, length(drive_t) - 1)
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

    m1 = (0.1 * (50.0 - Vs1) / (exp((50.0 - Vs1) / 10.0) - 1.0)) /
         ((0.1 * (50.0 - Vs1) / (exp((50.0 - Vs1) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs1) / 18.0))
    m2 = (0.1 * (50.0 - Vs2) / (exp((50.0 - Vs2) / 10.0) - 1.0)) /
         ((0.1 * (50.0 - Vs2) / (exp((50.0 - Vs2) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs2) / 18.0))
    m3 = (0.1 * (50.0 - Vs3) / (exp((50.0 - Vs3) / 10.0) - 1.0)) /
         ((0.1 * (50.0 - Vs3) / (exp((50.0 - Vs3) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs3) / 18.0))
    m4 = (0.1 * (50.0 - Vs4) / (exp((50.0 - Vs4) / 10.0) - 1.0)) /
         ((0.1 * (50.0 - Vs4) / (exp((50.0 - Vs4) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs4) / 18.0))

    du[1] = 0.0
    du[2] = 4.0 * m1^3 * h1 * (30.0 - V1) +
            0.3 * n1^4 * (-75.0 - V1) +
            0.01 * x1 * (30.0 - V1) +
            0.03 * Ca1 / (0.5 + Ca1) * (-75.0 - V1) +
            0.003 * (-40.0 - V1) -
            g41_eff * (V1 - 30.0) * s4 -
            params.g12 * (V1 + 80.0) * s21 / params.scale2 +
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
            params.g34 * (V3 + 80.0) * s43 / params.scale3 +
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

function simulate_publication_model_driven(params::ModelParams, control_gain::Float64, mode::ControlMode, drive_time_s::AbstractVector{<:Real}, drive_voltage_mv::AbstractVector{<:Real}; saveat_ms::Float64 = DRIVEN_TRACE_SAVEAT_MS)
    sim_params = mode_params(params, mode)
    drive_t = Float64.(drive_time_s)
    drive_v = Float64.(drive_voltage_mv)
    drive_t .-= first(drive_t)
    tspan_ms = (0.0, last(drive_t) * 1000.0)
    problem = ODEProblem(
        publication_driven_network_ode!,
        initial_state(),
        tspan_ms,
        (params = sim_params, control_gain = control_gain, mode = mode, drive_t = drive_t, drive_v = drive_v),
    )
    solution = solve(
        problem,
        RK4();
        saveat = saveat_ms,
        reltol = 1e-6,
        abstol = 1e-6,
        maxiters = 10^9,
    )
    time_s = solution.t ./ 1000.0
    states = Array(solution)
    driven_v0 = [interpolate_drive(t, drive_t, drive_v) for t in time_s]
    return (
        time_s = time_s,
        V0 = driven_v0,
        V1 = states[2, :],
        V2 = states[3, :],
        V3 = states[4, :],
        V4 = states[5, :],
    )
end

function representative_trace_initial_state(mode::ControlMode)
    if mode == presynaptic
        return copy(REPRESENTATIVE_TRACE_U0_PRESYNAPTIC)
    elseif mode == postsynaptic
        return copy(REPRESENTATIVE_TRACE_U0_POSTSYNAPTIC)
    end
    error("Unsupported representative-trace mode $(mode)")
end

function simulate_publication_model_driven_window(params::ModelParams, control_gain::Float64, mode::ControlMode, drive_time_s::AbstractVector{<:Real}, drive_voltage_mv::AbstractVector{<:Real}; saveat_ms::Float64 = DRIVEN_TRACE_SAVEAT_MS)
    drive_t = Float64.(drive_time_s)
    drive_v = Float64.(drive_voltage_mv)
    drive_t .-= first(drive_t)
    sim_params = mode_params(params, mode)
    u0 = representative_trace_initial_state(mode)
    problem = ODEProblem(
        publication_driven_network_ode!,
        u0,
        (0.0, last(drive_t) * 1000.0),
        (params = sim_params, control_gain = control_gain, mode = mode, drive_t = drive_t, drive_v = drive_v),
    )
    solution = solve(
        problem,
        RK4();
        saveat = saveat_ms,
        reltol = 1e-6,
        abstol = 1e-6,
        maxiters = 10^9,
    )
    time_s = solution.t ./ 1000.0
    states = Array(solution)
    driven_v0 = [interpolate_drive(t, drive_t, drive_v) for t in time_s]
    return (
        time_s = time_s,
        V0 = driven_v0,
        V1 = states[2, :],
        V2 = states[3, :],
        V3 = states[4, :],
        V4 = states[5, :],
    )
end

function simulate_compare_model_driven(mode_name::String, drive_time_s::AbstractVector{<:Real}, drive_voltage_mv::AbstractVector{<:Real}; step_ms::Float64 = 0.1, saveat_ms::Float64 = 2.0, alpha_scale::Float64 = 1.0, gm_scale::Float64 = 1.0, g_scale::Float64 = 1.0)
    cfg = compare_model_config(mode_name)
    drive_t = Float64.(drive_time_s)
    drive_v = Float64.(drive_voltage_mv)
    tmax_s = drive_t[end] - drive_t[1]
    drive_t .-= first(drive_t)

    alpha1 = 0.01
    beta1 = 0.002
    scale1 = alpha1 / (alpha1 + beta1)
    alpham_g = 0.005
    betam_g = 0.0001
    scalem_g = (alpham_g - betam_g) / alpham_g
    gm = 3.0
    alpham_a = 0.005
    betam_a = 0.0001
    scalem_a = (alpham_a - betam_a) / alpham_a
    gm_a = 3.0
    alphax = 0.011
    betax = 0.001
    alphai = 0.015
    betai = 0.0005
    scalei = (alphai - betai) / alphai
    alpha3 = 0.01
    beta3 = 0.002
    scale3 = alpha3 / (alpha3 + beta3)
    alpha2 = 0.01
    beta2 = 0.01
    scale2 = alpha2 / (alpha2 + beta2)
    gelec = 0.002
    g14 = 0.005
    g23 = 0.005
    g34 = 0.005
    g43 = 0.005
    g21 = 0.01
    g12 = 0.01
    Ca_shift4 = -25.0
    Ca_shift1 = -25.0
    Ca_shift3 = Ca_shift4
    Ca_shift2 = Ca_shift1
    x_shift = -4.0
    Vhh = -53.0

    V1 = -44.0
    V2 = -44.0
    V3 = -44.0
    V4 = -44.0
    x1 = 0.6
    x2 = 0.6
    x3 = 0.5
    x4 = 0.6
    Ca1 = 1.0
    Ca2 = 1.0
    Ca3 = 1.1
    Ca4 = 1.0
    h1 = 0.0
    h2 = 0.0
    h3 = 0.0
    h4 = 0.0
    n1 = 0.0
    n2 = 0.0
    n3 = 0.0
    n4 = 0.0
    y1 = 0.0
    y2 = 0.0
    s0 = 0.0
    s1 = 0.0
    s2 = 0.0
    s3 = 0.0
    s4 = 0.0
    s12 = 0.0
    s21 = 0.0
    s34 = 0.1
    s43 = 0.1
    sm_g = cfg.sm_g0
    sm_a = cfg.sm_a0

    nt = round(Int, tmax_s * 1000 / step_ms)
    isave = round(Int, saveat_ms / step_ms)
    nsave = round(Int, tmax_s * 1000 / saveat_ms) + 1
    drive_step = [interpolate_drive((i - 1) * step_ms / 1000.0, drive_t, drive_v) for i in 1:(nt + 1)]
    drive_save = [interpolate_drive((i - 1) * saveat_ms / 1000.0, drive_t, drive_v) for i in 1:nsave]

    time_s = zeros(Float64, nsave)
    vv0 = zeros(Float64, nsave)
    vv1 = zeros(Float64, nsave)
    save_idx = 1
    vv0[1] = drive_save[1]
    vv1[1] = V1

    for i in 1:nt
        V0 = drive_step[i]
        Vs0 = 127.0 * V0 / 105.0 + 8265.0 / 105.0
        Vs1 = 127.0 * V1 / 105.0 + 8265.0 / 105.0
        Vs2 = 127.0 * V2 / 105.0 + 8265.0 / 105.0
        Vs3 = 127.0 * V3 / 105.0 + 8265.0 / 105.0
        Vs4 = 127.0 * V4 / 105.0 + 8265.0 / 105.0

        m1 = (0.1 * (50.0 - Vs1) / (exp((50.0 - Vs1) / 10.0) - 1.0)) / ((0.1 * (50.0 - Vs1) / (exp((50.0 - Vs1) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs1) / 18.0))
        m2 = (0.1 * (50.0 - Vs2) / (exp((50.0 - Vs2) / 10.0) - 1.0)) / ((0.1 * (50.0 - Vs2) / (exp((50.0 - Vs2) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs2) / 18.0))
        m3 = (0.1 * (50.0 - Vs3) / (exp((50.0 - Vs3) / 10.0) - 1.0)) / ((0.1 * (50.0 - Vs3) / (exp((50.0 - Vs3) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs3) / 18.0))
        m4 = (0.1 * (50.0 - Vs4) / (exp((50.0 - Vs4) / 10.0) - 1.0)) / ((0.1 * (50.0 - Vs4) / (exp((50.0 - Vs4) / 10.0) - 1.0)) + 4.0 * exp((25.0 - Vs4) / 18.0))

        g41_eff = cfg.use_g_mod ? (cfg.g41 * g_scale * (1.0 + gm * gm_scale * sm_g / scalem_g)) : (cfg.g41 * g_scale)
        g32_eff = cfg.use_g_mod ? (cfg.g32 * g_scale * (1.0 + gm * gm_scale * sm_g / scalem_g)) : (cfg.g32 * g_scale)
        alpha_exc = cfg.use_a_mod ? (alphax * alpha_scale * (1.0 + gm_a * gm_scale * sm_a / scalem_a)) : (alphax * alpha_scale)

        V1_new = V1 + step_ms * (
            4.0 * m1^3 * h1 * (30.0 - V1) +
            0.3 * n1^4 * (-75.0 - V1) +
            0.01 * x1 * (30.0 - V1) +
            0.03 * Ca1 / (0.5 + Ca1) * (-75.0 - V1) +
            0.003 * (-40.0 - V1) -
            g41_eff * (V1 - 30.0) * s4 -
            g12 * (V1 + 80.0) * s21 / scale2 +
            gelec * (V4 - V1)
        )
        V2_new = V2 + step_ms * (
            4.0 * m2^3 * h2 * (30.0 - V2) +
            0.3 * n2^4 * (-75.0 - V2) +
            0.01 * x2 * (30.0 - V2) +
            0.03 * Ca2 / (0.5 + Ca2) * (-75.0 - V2) +
            0.003 * (-40.0 - V2) -
            g32_eff * (V2 - 30.0) * s3 -
            g12 * (V2 + 80.0) * s12 / scale2 +
            gelec * (V3 - V2)
        )
        V3_new = V3 + step_ms * (
            4.0 * m3^3 * h3 * (30.0 - V3) +
            0.3 * n3^4 * (-75.0 - V3) +
            0.01 * x3 * (30.0 - V3) +
            0.03 * Ca3 / (0.5 + Ca3) * (-75.0 - V3) +
            0.003 * (-40.0 - V3) -
            g23 * (V3 + 80.0) * s2 / scalei -
            g34 * (V3 + 80.0) * s43 / scale3 +
            gelec * (V2 - V3) -
            cfg.g0 * (V3 - 30.0) * s0 / scale1
        )
        V4_new = V4 + step_ms * (
            4.0 * m4^3 * h4 * (30.0 - V4) +
            0.3 * n4^4 * (-75.0 - V4) +
            0.01 * x4 * (30.0 - V4) +
            0.03 * Ca4 / (0.5 + Ca4) * (-75.0 - V4) +
            0.003 * (-40.0 - V4) -
            g14 * (V4 + 80.0) * s1 / scalei -
            g34 * (V4 + 80.0) * s34 / scale3 +
            gelec * (V1 - V4) -
            cfg.g0 * (V4 - 30.0) * s0 / scale1
        )

        Ca1_new = Ca1 + step_ms * (0.0003 * (0.0085 * x1 * (140.0 - V1 + Ca_shift1) - Ca1))
        Ca2_new = Ca2 + step_ms * (0.0003 * (0.0085 * x2 * (140.0 - V2 + Ca_shift2) - Ca2))
        Ca3_new = Ca3 + step_ms * (0.0003 * (0.0085 * x3 * (140.0 - V3 + Ca_shift3) - Ca3))
        Ca4_new = Ca4 + step_ms * (0.0003 * (0.0085 * x4 * (140.0 - V4 + Ca_shift4) - Ca4))

        x1_new = x1 + step_ms * (((1.0 / (exp(0.15 * (-V1 - 50.0 + x_shift)) + 1.0)) - x1) / 100.0)
        x2_new = x2 + step_ms * (((1.0 / (exp(0.15 * (-V2 - 50.0 + x_shift)) + 1.0)) - x2) / 100.0)
        x3_new = x3 + step_ms * (((1.0 / (exp(0.15 * (-V3 - 50.0 + x_shift)) + 1.0)) - x3) / 100.0)
        x4_new = x4 + step_ms * (((1.0 / (exp(0.15 * (-V4 - 50.0 + x_shift)) + 1.0)) - x4) / 100.0)

        h1_new = h1 + step_ms * (((1.0 - h1) * (0.07 * exp((25.0 - Vs1) / 20.0)) - h1 * (1.0 / (1.0 + exp((55.0 - Vs1) / 10.0)))) / 12.5)
        h2_new = h2 + step_ms * (((1.0 - h2) * (0.07 * exp((25.0 - Vs2) / 20.0)) - h2 * (1.0 / (1.0 + exp((55.0 - Vs2) / 10.0)))) / 12.5)
        h3_new = h3 + step_ms * (((1.0 - h3) * (0.07 * exp((25.0 - Vs3) / 20.0)) - h3 * (1.0 / (1.0 + exp((55.0 - Vs3) / 10.0)))) / 12.5)
        h4_new = h4 + step_ms * (((1.0 - h4) * (0.07 * exp((25.0 - Vs4) / 20.0)) - h4 * (1.0 / (1.0 + exp((55.0 - Vs4) / 10.0)))) / 12.5)

        n1_new = n1 + step_ms * (((1.0 - n1) * (0.01 * (55.0 - Vs1) / (exp((55.0 - Vs1) / 10.0) - 1.0)) - n1 * (0.125 * exp((45.0 - Vs1) / 80.0))) / 12.5)
        n2_new = n2 + step_ms * (((1.0 - n2) * (0.01 * (55.0 - Vs2) / (exp((55.0 - Vs2) / 10.0) - 1.0)) - n2 * (0.125 * exp((45.0 - Vs2) / 80.0))) / 12.5)
        n3_new = n3 + step_ms * (((1.0 - n3) * (0.01 * (55.0 - Vs3) / (exp((55.0 - Vs3) / 10.0) - 1.0)) - n3 * (0.125 * exp((45.0 - Vs3) / 80.0))) / 12.5)
        n4_new = n4 + step_ms * (((1.0 - n4) * (0.01 * (55.0 - Vs4) / (exp((55.0 - Vs4) / 10.0) - 1.0)) - n4 * (0.125 * exp((45.0 - Vs4) / 80.0))) / 12.5)

        y1_new = y1 + step_ms * (0.5 * ((1.0 / (1.0 + exp(10.0 * (V1 - Vhh)))) - y1) / (7.1 + 10.4 / (1.0 + exp((V1 + 68.0) / 2.2))))
        y2_new = y2 + step_ms * (0.5 * ((1.0 / (1.0 + exp(10.0 * (V2 - Vhh)))) - y2) / (7.1 + 10.4 / (1.0 + exp((V2 + 68.0) / 2.2))))

        act0 = 1.0 / (1.0 + exp(-20.0 * (V0 + 20.0)))
        act1 = 1.0 / (1.0 + exp(-20.0 * (V1 + 20.0)))
        act2 = 1.0 / (1.0 + exp(-20.0 * (V2 + 20.0)))
        act3 = 1.0 / (1.0 + exp(-20.0 * (V3 + 20.0)))
        act4 = 1.0 / (1.0 + exp(-20.0 * (V4 + 20.0)))

        s0_new = s0 + step_ms * (alpha1 * (1.0 - s0) * act0 - beta1 * s0)
        sm_g_new = cfg.use_g_mod ? sm_g + step_ms * (alpham_g * sm_g * (1.0 - sm_g) * act0 - betam_g * (sm_g - 0.0001)) : sm_g
        sm_a_new = cfg.use_a_mod ? sm_a + step_ms * (alpham_a * sm_a * (1.0 - sm_a) * act0 - betam_a * (sm_a - 0.0001)) : sm_a
        s1_new = s1 + step_ms * (alphai * s1 * (1.0 - s1) * act1 - betai * (s1 - 0.0001))
        s2_new = s2 + step_ms * (alphai * s2 * (1.0 - s2) * act2 - betai * (s2 - 0.0001))
        s3_new = s3 + step_ms * (alpha_exc * s3 * (1.0 - s3) * act3 - betax * (s3 - 0.0001))
        s4_new = s4 + step_ms * (alpha_exc * s4 * (1.0 - s4) * act4 - betax * (s4 - 0.0001))
        s12_new = s12 + step_ms * (alpha2 * (1.0 - s12) * act1 - beta2 * s12)
        s21_new = s21 + step_ms * (alpha2 * (1.0 - s21) * act2 - beta2 * s21)
        s34_new = s34 + step_ms * (alpha3 * (1.0 - s34) * act3 - beta3 * s34)
        s43_new = s43 + step_ms * (alpha3 * (1.0 - s43) * act4 - beta3 * s43)

        V1, V2, V3, V4 = V1_new, V2_new, V3_new, V4_new
        x1, x2, x3, x4 = x1_new, x2_new, x3_new, x4_new
        Ca1, Ca2, Ca3, Ca4 = Ca1_new, Ca2_new, Ca3_new, Ca4_new
        h1, h2, h3, h4 = h1_new, h2_new, h3_new, h4_new
        n1, n2, n3, n4 = n1_new, n2_new, n3_new, n4_new
        y1, y2 = y1_new, y2_new
        s0, s1, s2, s3, s4 = s0_new, s1_new, s2_new, s3_new, s4_new
        s12, s21, s34, s43 = s12_new, s21_new, s34_new, s43_new
        sm_g, sm_a = sm_g_new, sm_a_new

        if i % isave == 0
            save_idx += 1
            time_s[save_idx] = i * step_ms / 1000.0
            vv0[save_idx] = drive_save[save_idx]
            vv1[save_idx] = V1
        end
    end

    return DataFrame(
        mode = fill(mode_name, save_idx),
        time_s = time_s[1:save_idx],
        V0 = vv0[1:save_idx],
        V1 = vv1[1:save_idx],
    )
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

    spike_times = upward_spike_times(time_s, v0; threshold_mv = -20.0, refractory_s = 0.05)
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

function ensure_master_model_cache()
    if isfile(MASTER_MODEL_POINTS_CSV) && isfile(MASTER_MODEL_TRACES_CSV)
        return CSV.read(MASTER_MODEL_TRACES_CSV, DataFrame), CSV.read(MASTER_MODEL_POINTS_CSV, DataFrame)
    end

    pre_trace = simulate_compare_model("presynaptic")
    post_trace = simulate_compare_model("postsynaptic")
    traces = vcat(pre_trace, post_trace)
    points = vcat(matched_si1_vs_si2_points(pre_trace), matched_si1_vs_si2_points(post_trace))
    CSV.write(MASTER_MODEL_TRACES_CSV, traces)
    CSV.write(MASTER_MODEL_POINTS_CSV, points)
    return traces, points
end

function ensure_master_driven_trace_cache(drive_time_s::AbstractVector{<:Real}, drive_voltage_mv::AbstractVector{<:Real}; force::Bool = false)
    if !force && isfile(MASTER_DRIVEN_TRACES_CSV)
        cached = CSV.read(MASTER_DRIVEN_TRACES_CSV, DataFrame)
        if finite_trace_cache(cached)
            return cached
        end
    end
    summary = publication_summary_subset()
    params = publication_params()
    pre_gain = representative_gain(summary, "presynaptic")
    post_gain = representative_gain(summary, "postsynaptic")
    pre = simulation_dataframe(simulate_publication_model_driven_window(params, pre_gain, presynaptic, drive_time_s, drive_voltage_mv), "presynaptic")
    post = simulation_dataframe(simulate_publication_model_driven_window(params, post_gain, postsynaptic, drive_time_s, drive_voltage_mv), "postsynaptic")
    traces = vcat(pre, post)
    finite_trace_cache(traces) || error("Driven representative trace cache contains non-finite values after regeneration.")
    CSV.write(MASTER_DRIVEN_TRACES_CSV, traces)
    return traces
end

function prepare_master_figure_data(; force::Bool = false)
    mkpath(OUTPUT_DIR)

    if force || !isfile(PUB_SCAN_SUMMARY_CSV) || !isfile(PUB_SCAN_RAW_CSV)
        raw, summary = generate_publication_data()
        save_publication_data(raw, summary)
    end

    if force || !isfile(MASTER_MODEL_POINTS_CSV) || !isfile(MASTER_MODEL_TRACES_CSV)
        pre_trace = simulate_compare_model("presynaptic")
        post_trace = simulate_compare_model("postsynaptic")
        traces = vcat(pre_trace, post_trace)
        points = vcat(matched_si1_vs_si2_points(pre_trace), matched_si1_vs_si2_points(post_trace))
        CSV.write(MASTER_MODEL_TRACES_CSV, traces)
        CSV.write(MASTER_MODEL_POINTS_CSV, points)
    end

    bio_t1, bio_v1 = load_two_column_txt(BIO_SI1_TXT)
    bio_t1 = bio_t1 .- first(bio_t1)
    ensure_master_driven_trace_cache(bio_t1, bio_v1; force = force)
end

function load_master_figure_data()
    missing = String[]
    isfile(PUB_SCAN_SUMMARY_CSV) || push!(missing, PUB_SCAN_SUMMARY_CSV)
    isfile(PUB_SCAN_RAW_CSV) || push!(missing, PUB_SCAN_RAW_CSV)
    isfile(MASTER_MODEL_POINTS_CSV) || push!(missing, MASTER_MODEL_POINTS_CSV)
    isfile(MASTER_MODEL_TRACES_CSV) || push!(missing, MASTER_MODEL_TRACES_CSV)
    isfile(MASTER_DRIVEN_TRACES_CSV) || push!(missing, MASTER_DRIVEN_TRACES_CSV)
    isempty(missing) || error("Missing master-figure cache files. Run prepare_master_figure_data() or `julia --project=jack_figs .\\jack_figs\\jl\\neuromod_master_figure.jl prepare` first.\n" * join(missing, "\n"))

    return (
        traces = CSV.read(MASTER_MODEL_TRACES_CSV, DataFrame),
        model_points = CSV.read(MASTER_MODEL_POINTS_CSV, DataFrame),
        summary = publication_summary_subset(),
        raw = publication_raw_subset(),
        driven_traces = CSV.read(MASTER_DRIVEN_TRACES_CSV, DataFrame),
    )
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
    return
    limits!(ax, 0, 10, 0, 10)

    node_xy = Dict(
        "Si2L" => Point2f(2, 8),
        "Si2R" => Point2f(8, 8),
        "Si3L" => Point2f(2, 4.6),
        "Si3R" => Point2f(8, 4.6),
        "Si1" => Point2f(5, 1.3),
    )
    node_color = Dict(
        "Si2L" => RGBf(0.15, 0.35, 0.90),
        "Si2R" => RGBf(0.15, 0.35, 0.90),
        "Si3L" => RGBf(0.88, 0.22, 0.20),
        "Si3R" => RGBf(0.88, 0.22, 0.20),
        "Si1" => RGBf(0.35, 0.35, 0.35),
    )

    function connect(a, b; color=:black, width=3, style=nothing)
        p1 = node_xy[a]
        p2 = node_xy[b]
        GLMakie.lines!(ax, [p1[1], p2[1]], [p1[2], p2[2]]; color = color, linewidth = width, linestyle = style)
    end

    connect("Si2L", "Si2R"; color = (:black, 0.45), width = 2)
    connect("Si3L", "Si3R"; color = (:black, 0.45), width = 2)
    connect("Si2L", "Si3R"; color = MODE_COLORS["postsynaptic"], width = 4)
    connect("Si2R", "Si3L"; color = MODE_COLORS["postsynaptic"], width = 4)
    connect("Si3L", "Si2L"; color = MODE_COLORS["presynaptic"], width = 4)
    connect("Si3R", "Si2R"; color = MODE_COLORS["presynaptic"], width = 4)
    connect("Si1", "Si3L"; color = (:black, 0.8), width = 3)
    connect("Si1", "Si3R"; color = (:black, 0.8), width = 3)

    GLMakie.lines!(ax, [5.0, 6.8], [1.3, 6.2]; color = (:black, 0.8), linewidth = 2, linestyle = :dash)
    GLMakie.text!(ax, 7.0, 6.5; text = "modulates\nSi3→Si2", align = (:left, :bottom), fontsize = 16)
    GLMakie.text!(ax, 7.0, 5.6; text = "α or g", align = (:left, :bottom), fontsize = 16, color = MODE_COLORS["presynaptic"])

    for (name, pt) in node_xy
        GLMakie.scatter!(ax, [pt[1]], [pt[2]]; color = node_color[name], markersize = 45)
        GLMakie.text!(ax, pt[1], pt[2]; text = name, align = (:center, :center), fontsize = 18, color = :white)
    end
    ax.title = "Circuit"
end

function plot_stacked_trace_panel!(ax::Axis, upper_time_s::AbstractVector{<:Real}, upper_v::AbstractVector{<:Real}, lower_time_s::AbstractVector{<:Real}, lower_v::AbstractVector{<:Real}, upper_label::String, lower_label::String; upper_color = RGBf(0.35, 0.35, 0.35), lower_color = RGBf(0.15, 0.35, 0.90), title::String, time_bar_s::Float64 = 20.0, voltage_bar_mv::Float64 = 50.0)
    offset = 95.0
    upper_t = Float64.(upper_time_s)
    lower_t = Float64.(lower_time_s)
    upper_vals = Float64.(upper_v)
    lower_vals = Float64.(lower_v)
    GLMakie.lines!(ax, upper_t, upper_vals .+ offset; color = upper_color, linewidth = 1.5)
    GLMakie.lines!(ax, lower_t, lower_vals; color = lower_color, linewidth = 1.5)
    x_min = min(first(upper_t), first(lower_t))
    x_max = max(last(upper_t), last(lower_t))
    y_min = min(minimum(lower_vals), minimum(upper_vals .+ offset))
    y_max = max(maximum(lower_vals), maximum(upper_vals .+ offset))
    GLMakie.text!(ax, x_min + 0.015 * (x_max - x_min), offset + 38; text = upper_label, align = (:left, :center), fontsize = 17, color = upper_color)
    GLMakie.text!(ax, x_min + 0.015 * (x_max - x_min), 38; text = lower_label, align = (:left, :center), fontsize = 17, color = lower_color)
    ax.title = title
    ax.xlabel = ""
    ax.ylabel = "V (mV)"
    ax.xticksvisible = false
    ax.yticksvisible = false
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    hidespines!(ax, :t, :r)
    add_l_scale_bar!(ax, x_min, x_max, y_min, y_max; time_len_s = time_bar_s, volt_len_mv = voltage_bar_mv)
end

function plot_three_trace_panel!(ax::Axis, time_s::AbstractVector{<:Real}, top_v::AbstractVector{<:Real}, mid_v::AbstractVector{<:Real}, bot_v::AbstractVector{<:Real}, top_label::String, mid_label::String, bot_label::String; top_color = RGBf(0.35, 0.35, 0.35), mid_color = MODE_COLORS["presynaptic"], bot_color = MODE_COLORS["postsynaptic"], title::String, time_bar_s::Float64 = 20.0, voltage_bar_mv::Float64 = 50.0)
    offset_top = 190.0
    offset_mid = 95.0
    t = Float64.(time_s)
    top_vals = Float64.(top_v)
    mid_vals = Float64.(mid_v)
    bot_vals = Float64.(bot_v)
    GLMakie.lines!(ax, t, top_vals .+ offset_top; color = top_color, linewidth = 1.5)
    GLMakie.lines!(ax, t, mid_vals .+ offset_mid; color = mid_color, linewidth = 1.5)
    GLMakie.lines!(ax, t, bot_vals; color = bot_color, linewidth = 1.5)
    x_min = first(t)
    x_max = last(t)
    y_min = minimum(bot_vals)
    y_max = maximum(top_vals .+ offset_top)
    GLMakie.text!(ax, x_min + 0.015 * (x_max - x_min), offset_top + 38; text = top_label, align = (:left, :center), fontsize = 17, color = top_color)
    GLMakie.text!(ax, x_min + 0.015 * (x_max - x_min), offset_mid + 38; text = mid_label, align = (:left, :center), fontsize = 17, color = mid_color)
    GLMakie.text!(ax, x_min + 0.015 * (x_max - x_min), 38; text = bot_label, align = (:left, :center), fontsize = 17, color = bot_color)
    ax.title = title
    ax.xlabel = ""
    ax.ylabel = "V (mV)"
    ax.xticksvisible = false
    ax.yticksvisible = false
    ax.xticklabelsvisible = false
    ax.yticklabelsvisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    hidespines!(ax, :t, :r)
    add_l_scale_bar!(ax, x_min, x_max, y_min, y_max; time_len_s = time_bar_s, volt_len_mv = voltage_bar_mv)
end

function plot_four_trace_panel!(ax::Axis, t1::AbstractVector{<:Real}, v1::AbstractVector{<:Real}, t2::AbstractVector{<:Real}, v2::AbstractVector{<:Real}, t3::AbstractVector{<:Real}, v3::AbstractVector{<:Real}, t4::AbstractVector{<:Real}, v4::AbstractVector{<:Real}, label1::String, label2::String, label3::String, label4::String; color1 = RGBf(0.35, 0.35, 0.35), color2 = RGBf(0.1, 0.1, 0.1), color3 = MODE_COLORS["presynaptic"], color4 = MODE_COLORS["postsynaptic"], title::String, time_bar_s::Float64 = 20.0, voltage_bar_mv::Float64 = 50.0)
    offsets = [285.0, 190.0, 95.0, 0.0]
    time1 = Float64.(t1)
    time2 = Float64.(t2)
    time3 = Float64.(t3)
    time4 = Float64.(t4)
    vals1 = Float64.(v1)
    vals2 = Float64.(v2)
    vals3 = Float64.(v3)
    vals4 = Float64.(v4)
    GLMakie.lines!(ax, time1, vals1 .+ offsets[1]; color = color1, linewidth = 1.4)
    GLMakie.lines!(ax, time2, vals2 .+ offsets[2]; color = color2, linewidth = 1.4)
    GLMakie.lines!(ax, time3, vals3 .+ offsets[3]; color = color3, linewidth = 1.4)
    GLMakie.lines!(ax, time4, vals4 .+ offsets[4]; color = color4, linewidth = 1.4)
    x_min = minimum((first(time1), first(time2), first(time3), first(time4)))
    x_max = maximum((last(time1), last(time2), last(time3), last(time4)))
    y_min = minimum(vals4)
    y_max = maximum(vals1 .+ offsets[1])
    label_x = x_min + 0.015 * (x_max - x_min)
    GLMakie.text!(ax, label_x, offsets[1] + 38; text = label1, align = (:left, :center), fontsize = 16, color = color1)
    GLMakie.text!(ax, label_x, offsets[2] + 38; text = label2, align = (:left, :center), fontsize = 16, color = color2)
    GLMakie.text!(ax, label_x, offsets[3] + 38; text = label3, align = (:left, :center), fontsize = 16, color = color3)
    GLMakie.text!(ax, label_x, offsets[4] + 38; text = label4, align = (:left, :center), fontsize = 16, color = color4)
    ax.title = title
    ax.xlabel = ""
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
    add_l_scale_bar!(ax, x_min, x_max, y_min, y_max; time_len_s = time_bar_s, volt_len_mv = voltage_bar_mv)
end

function plot_gain_panel!(ax::Axis, summary::DataFrame, raw::DataFrame, phase::String)
    ylabel = phase == "pre" ? "Freq on (Hz)" : "Freq off (Hz)"
    title = phase == "pre" ? "During Drive" : "After Drive Removal"
    ax.xlabel = "Gain"
    ax.ylabel = ylabel
    ax.title = title
    GLMakie.xlims!(ax, 0, 6)
    ylo = phase == "post" ? -0.005 : 0.0
    GLMakie.ylims!(ax, ylo, 0.21)

    legend_entries = Any[]
    legend_labels = String[]
    for mode_name in ["presynaptic", "postsynaptic"]
        color = MODE_COLORS[mode_name]
        marker = MODE_MARKERS[mode_name]
        sub_summary = sort(summary[summary[!, :mode] .== mode_name, :], :control_gain)
        sub_raw = raw[(raw[!, :mode] .== mode_name) .& (raw[!, :phase] .== phase), :]
        bgx, bgy = finite_xy(sub_raw[!, :control_gain], sub_raw[!, :frequency_hz])
        xcol = phase == "pre" ? :mean_pre_frequency_hz : :mean_post_frequency_hz
        x, y = finite_xy(sub_summary[!, :control_gain], sub_summary[!, xcol])
        if !isempty(bgx)
            GLMakie.scatter!(ax, bgx, bgy; color = (color, 0.24), marker = :circle, markersize = 7, strokewidth = 0)
        end
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
        sub = model_points[model_points[!, :mode] .== mode_name, :]
        handle = GLMakie.scatter!(
            ax,
            sub[!, :si1_frequency_hz],
            sub[!, :burst_frequency_hz];
            color = (MODE_COLORS[mode_name], 0.84),
            marker = MODE_MARKERS[mode_name],
            markersize = 11,
            strokecolor = :black,
            strokewidth = 0.8,
        )
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

function build_master_figure(; ensure_cache::Bool = false)
    ensure_cache && prepare_master_figure_data()
    caches = load_master_figure_data()
    traces = caches.traces
    model_points = caches.model_points
    summary = caches.summary
    raw = caches.raw

    bio_t1, bio_v1 = load_two_column_txt(BIO_SI1_TXT)
    bio_t2, bio_v2 = load_two_column_txt(BIO_SI2_TXT)
    bio_t1 = bio_t1 .- first(bio_t1)
    bio_t2 = bio_t2 .- first(bio_t2)
    driven_traces = caches.driven_traces
    pre_driven = driven_traces[driven_traces[!, :mode] .== "presynaptic", :]
    post_driven = driven_traces[driven_traces[!, :mode] .== "postsynaptic", :]
    model_t = Float64.(pre_driven[!, :time_s])
    model_pre_v1 = Float64.(pre_driven[!, :V1])
    model_post_v1 = Float64.(post_driven[!, :V1])

    trace_start = 0.0
    trace_stop = min(maximum(bio_t1), maximum(model_t))
    bio1_keep = (bio_t1 .>= trace_start) .& (bio_t1 .<= trace_stop)
    bio2_keep = (bio_t2 .>= trace_start) .& (bio_t2 .<= trace_stop)
    model_keep = (model_t .>= trace_start) .& (model_t .<= trace_stop)

    bio_t1_view = bio_t1[bio1_keep] .- trace_start
    bio_v1_view = bio_v1[bio1_keep]
    bio_t2_view = bio_t2[bio2_keep] .- trace_start
    bio_v2_view = bio_v2[bio2_keep]
    model_t_view = model_t[model_keep] .- trace_start
    model_pre_v1_view = model_pre_v1[model_keep]
    model_post_v1_view = model_post_v1[model_keep]

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
    plot_four_trace_panel!(axE, bio_t1_view, bio_v1_view, bio_t2_view, bio_v2_view, model_t_view, model_pre_v1_view, model_t_view, model_post_v1_view, "Si1 biology", "Si2 biology", "Si2 presynaptic", "Si2 postsynaptic"; title = "Representative Traces", time_bar_s = 20.0)

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

    CairoMakie.save(MASTER_PNG, fig)
    CairoMakie.save(PAPER_FIG5_JPG, fig)
    return fig
end

const MASTER_NOAUTORUN = get(ENV, "NEUROMOD_MASTER_NOAUTORUN", "0") == "1"

function main(args = ARGS)
    mode = isempty(args) ? "both" : lowercase(args[1])
    if mode == "prepare"
        prepare_master_figure_data()
        println("Prepared master-figure caches:")
        println("  $(PUB_SCAN_RAW_CSV)")
        println("  $(PUB_SCAN_SUMMARY_CSV)")
        println("  $(MASTER_MODEL_TRACES_CSV)")
        println("  $(MASTER_MODEL_POINTS_CSV)")
        println("  $(MASTER_DRIVEN_TRACES_CSV)")
        return
    elseif mode == "plot"
        fig = build_master_figure(; ensure_cache = false)
        isinteractive() && display(fig)
        println("Saved $(MASTER_PNG)")
        return
    elseif mode == "both"
        prepare_master_figure_data()
        fig = build_master_figure(; ensure_cache = false)
        isinteractive() && display(fig)
        println("Saved $(MASTER_PNG)")
        return
    end
    error("Unknown mode `$(mode)`. Use `prepare`, `plot`, or `both`.")
end

if !MASTER_NOAUTORUN
    main()
end
