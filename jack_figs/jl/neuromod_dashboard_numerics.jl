using CSV
using DataFrames
using OrdinaryDiffEq
using Random
using Statistics

include(joinpath(@__DIR__, "bio_driven_traces.jl"))

const DASHBOARD_OUTPUT_DIR = joinpath(@__DIR__, "outputs")
const DASHBOARD_POINTS_CSV = joinpath(DASHBOARD_OUTPUT_DIR, "dashboard_scan_points.csv")
const DASHBOARD_SUMMARY_CSV = joinpath(DASHBOARD_OUTPUT_DIR, "dashboard_scan_summary.csv")
const DASHBOARD_PARAM_DEFAULTS_CSV = joinpath(DASHBOARD_OUTPUT_DIR, "dashboard_param_defaults.csv")
const DASHBOARD_CLEAN_PNG = joinpath(DASHBOARD_OUTPUT_DIR, "neuromod_dashboard_clean.png")
const DASHBOARD_CIRCUIT_SOURCE = joinpath(dirname(@__DIR__), "s1_Fig2A_drives_cpg_als_circuit.jpg")
const DASHBOARD_CIRCUIT_CROP = joinpath(DASHBOARD_OUTPUT_DIR, "legacy_circuit_crop2.jpg")
const BIO_FIG2C_TSV = joinpath(dirname(@__DIR__), "akira", "xlsx_cache", "Fig_2C.tsv")
const REPRESENTATIVE_BIOLOGY_TIME_SCALE = 2.0

const ORDER_LABELS = ["left-to-right", "right-to-left", "random"]
const MODE_LABELS = ["presynaptic", "postsynaptic"]

function load_dashboard_raw()
    isfile(DASHBOARD_POINTS_CSV) && return CSV.read(DASHBOARD_POINTS_CSV, DataFrame)
    isfile(PRESYNAPTIC_POINTS_CSV) && isfile(POSTSYNAPTIC_POINTS_CSV) &&
        return vcat(CSV.read(PRESYNAPTIC_POINTS_CSV, DataFrame), CSV.read(POSTSYNAPTIC_POINTS_CSV, DataFrame))
    return empty_publication_raw()
end

function load_dashboard_summary()
    isfile(DASHBOARD_SUMMARY_CSV) && return CSV.read(DASHBOARD_SUMMARY_CSV, DataFrame)
    isfile(PRESYNAPTIC_SUMMARY_CSV) && isfile(POSTSYNAPTIC_SUMMARY_CSV) &&
        return vcat(CSV.read(PRESYNAPTIC_SUMMARY_CSV, DataFrame), CSV.read(POSTSYNAPTIC_SUMMARY_CSV, DataFrame))
    return empty_publication_summary()
end

function ordered_gain_grid(gmin::Float64, gmax::Float64, n::Int, order::String)
    lo, hi = min(gmin, gmax), max(gmin, gmax)
    values = collect(range(lo, hi; length = max(n, 2)))
    if order == "right-to-left"
        reverse!(values)
    elseif order == "random"
        shuffle!(MersenneTwister(), values)
    elseif order != "left-to-right"
        error("Unknown scan order $(order)")
    end
    return Float64.(values)
end

function dashboard_config(max_time_s::Float64, saveat_ms::Float64; spike_threshold_mv::Float64, spike_refractory_s::Float64, burst_factor::Float64)
    base = calibrated_scan_config()
    SweepConfig(
        max_time_s,
        saveat_ms,
        spike_threshold_mv,
        spike_refractory_s,
        burst_factor,
        base.onset_timeout_s,
        base.transient_bursts,
        base.measured_cycles,
    )
end

function dashboard_params(mode::ControlMode; x_shift_si2::Float64, x_shift_si3::Float64, presyn_g0::Float64, postsyn_g0::Float64,
    alpha1::Float64, beta1::Float64, alpham::Float64, betam::Float64,
    ca_shift_si2::Float64, ca_shift_si3::Float64,
    presynaptic_base_g::Float64, direct_post_base_g::Float64, presyn_alphax::Float64, postsyn_alphax::Float64,
    si2_mutual_inhib_g::Float64, si2_mutual_inhib_alpha::Float64, si2_mutual_inhib_beta::Float64,
    si3_mutual_inhib_g::Float64, si3_mutual_inhib_alpha::Float64, si3_mutual_inhib_beta::Float64, t1_ms::Float64)
    params = calibrated_params(mode)
    return updated_params(
        params;
        x_shift = x_shift_si2,
        x_shift1 = x_shift_si2,
        x_shift2 = x_shift_si2,
        x_shift3 = x_shift_si3,
        x_shift4 = x_shift_si3,
        g0 = mode == presynaptic ? presyn_g0 : postsyn_g0,
        alpha1 = alpha1,
        beta1 = beta1,
        scale1 = alpha1 / (alpha1 + beta1),
        alpham = alpham,
        betam = betam,
        scale_sm = (alpham - betam) / alpham,
        presynaptic_base_alpha = presyn_alphax,
        Ca_shift1 = ca_shift_si2,
        Ca_shift2 = ca_shift_si2,
        Ca_shift3 = ca_shift_si3,
        Ca_shift4 = ca_shift_si3,
        presynaptic_base_g = presynaptic_base_g,
        direct_post_base_g = direct_post_base_g,
        alphax = postsyn_alphax,
        g12 = si2_mutual_inhib_g,
        g21 = si2_mutual_inhib_g,
        alpha2 = si2_mutual_inhib_alpha,
        beta2 = si2_mutual_inhib_beta,
        scale2 = si2_mutual_inhib_alpha / (si2_mutual_inhib_alpha + si2_mutual_inhib_beta),
        g34 = si3_mutual_inhib_g,
        g43 = si3_mutual_inhib_g,
        alpha3 = si3_mutual_inhib_alpha,
        beta3 = si3_mutual_inhib_beta,
        scale3 = si3_mutual_inhib_alpha / (si3_mutual_inhib_alpha + si3_mutual_inhib_beta),
        t1_ms = t1_ms,
    )
end

function row_for_selection(summary::DataFrame, mode_name::String, gain::Float64)
    sub = summary[summary.mode .== mode_name, :]
    isempty(sub) && return nothing
    idx = argmin(abs.(Float64.(sub.control_gain) .- gain))
    return sub[idx, :]
end

function run_scan_point(params::ModelParams, mode::ControlMode, gain::Float64, config::SweepConfig)
    sim = run_simulation(params, gain, mode, config)
    raw = publication_raw_points(sim, gain, mode, params, config)
    summary = publication_summary_row(sim, raw, gain, mode)
    summary[!, :si1_excitatory_g] .= params.g0
    summary[!, :alpha1] .= params.alpha1
    summary[!, :beta1] .= params.beta1
    summary[!, :alpham] .= params.alpham
    summary[!, :betam] .= params.betam
    summary[!, :ca_shift_si2] .= params.Ca_shift1
    summary[!, :ca_shift_si3] .= params.Ca_shift3
    summary[!, :x_shift_si2] .= params.x_shift1
    summary[!, :x_shift_si3] .= params.x_shift3
    summary[!, :presynaptic_base_g] .= params.presynaptic_base_g
    summary[!, :direct_post_base_g] .= params.direct_post_base_g
    summary[!, :presynaptic_base_alpha] .= params.presynaptic_base_alpha
    summary[!, :postsynaptic_base_alpha] .= params.alphax
    summary[!, :alphax] .= params.alphax
    summary[!, :si2_mutual_inhib_g] .= params.g12
    summary[!, :si2_mutual_inhib_alpha] .= params.alpha2
    summary[!, :si2_mutual_inhib_beta] .= params.beta2
    summary[!, :si3_mutual_inhib_g] .= params.g34
    summary[!, :si3_mutual_inhib_alpha] .= params.alpha3
    summary[!, :si3_mutual_inhib_beta] .= params.beta3
    summary[!, :t1_ms] .= params.t1_ms
    return raw, summary
end

function load_biological_fig2c(; time_scale::Float64 = REPRESENTATIVE_BIOLOGY_TIME_SCALE)
    df = DataFrame(CSV.File(BIO_FIG2C_TSV; delim = '\t', normalizenames = false))
    x = Float64[]
    y = Float64[]
    frequency_scale = 1.0 / time_scale
    for row in eachrow(df)
        xv = row[1]
        yv = row[2]
        if ismissing(xv) || ismissing(yv) || xv == "" || yv == ""
            continue
        end
        push!(x, (xv isa AbstractString ? parse(Float64, xv) : Float64(xv)) * frequency_scale)
        push!(y, (yv isa AbstractString ? parse(Float64, yv) : Float64(yv)) * frequency_scale)
    end
    return x, y
end

function load_biological_trace_points(; time_scale::Float64 = REPRESENTATIVE_BIOLOGY_TIME_SCALE, spike_threshold_mv::Float64 = -20.0, spike_refractory_s::Float64 = 0.02, burst_factor::Float64 = 2.0)
    si1_t, si1_v, si2_t, si2_v = load_biology_pair(; time_scale = time_scale)
    t_start = max(first(si1_t), first(si2_t))
    t_stop = min(last(si1_t), last(si2_t))
    keep = (si2_t .>= t_start) .& (si2_t .<= t_stop)
    common_t = Float64.(si2_t[keep])
    isempty(common_t) && return DataFrame(mode = String[], si1_frequency_hz = Float64[], burst_frequency_hz = Float64[])
    trace = DataFrame(
        mode = fill("biology", length(common_t)),
        time_s = common_t,
        V0 = [interpolate_drive(t, si1_t, si1_v) for t in common_t],
        V1 = Float64.(si2_v[keep]),
    )
    return matched_si1_vs_si2_points(trace; spike_threshold_mv = spike_threshold_mv, spike_refractory_s = spike_refractory_s, burst_factor = burst_factor)
end

function representative_model_points(traces::Dict{String, DataFrame}; spike_threshold_mv::Float64 = -20.0, spike_refractory_s::Float64 = 0.02, burst_factor::Float64 = 2.0)
    parts = DataFrame[]
    for mode_name in MODE_LABELS
        tr = traces[mode_name]
        isempty(tr) && continue
        push!(parts, matched_si1_vs_si2_points(tr; spike_threshold_mv = spike_threshold_mv, spike_refractory_s = spike_refractory_s, burst_factor = burst_factor))
    end
    isempty(parts) && return DataFrame(mode = String[], si1_frequency_hz = Float64[], burst_frequency_hz = Float64[])
    return vcat(parts...)
end

function empty_trace_df()
    DataFrame(mode = String[], control_gain = Float64[], time_s = Float64[], V0 = Float64[], V1 = Float64[], V3 = Float64[], V4 = Float64[])
end

function empty_representative_traces()
    Dict("presynaptic" => empty_trace_df(), "postsynaptic" => empty_trace_df())
end

function load_representative_traces()
    traces = empty_representative_traces()
    isfile(BIO_DRIVEN_TRACES_CSV) || return traces
    df = CSV.read(BIO_DRIVEN_TRACES_CSV, DataFrame)
    for mode_name in MODE_LABELS
        sub = df[df.mode .== mode_name, :]
        if !isempty(sub)
            traces[mode_name] = DataFrame(
                mode = String.(sub.mode),
                control_gain = fill(mode_name == "presynaptic" ? PRESYNAPTIC_ANCHOR_GAIN : POSTSYNAPTIC_ANCHOR_GAIN, nrow(sub)),
                time_s = Float64.(sub.time_s) .* REPRESENTATIVE_BIOLOGY_TIME_SCALE,
                V0 = Float64.(sub.V0),
                V1 = Float64.(sub.V1),
                V3 = fill(NaN, nrow(sub)),
                V4 = fill(NaN, nrow(sub)),
            )
        end
    end
    return traces
end

function load_biology_pair(; time_scale::Float64 = REPRESENTATIVE_BIOLOGY_TIME_SCALE)
    si1_t, si1_v = load_two_column_txt(BIO_SI1_TXT)
    si2_t, si2_v = load_two_column_txt(BIO_SI2_TXT)
    si1_t .-= first(si1_t)
    si2_t .-= first(si2_t)
    si1_t .*= time_scale
    si2_t .*= time_scale
    return si1_t, si1_v, si2_t, si2_v
end

function apply_calcium_initial_conditions(u0::Vector{Float64}, calcium_ics)
    isnothing(calcium_ics) && return u0
    out = copy(u0)
    out[11] = Float64(calcium_ics.Ca0)
    out[12] = Float64(calcium_ics.Ca1)
    out[13] = Float64(calcium_ics.Ca2)
    out[14] = Float64(calcium_ics.Ca3)
    out[15] = Float64(calcium_ics.Ca4)
    return out
end

function simulate_driven_state_trace(params::ModelParams, gain::Float64, mode::ControlMode,
    drive_time_s::Vector{Float64}, drive_voltage_mv::Vector{Float64}; saveat_ms::Float64, display_time_scale::Float64 = 1.0, calcium_ics = nothing)
    sim_params = mode_params(params, mode)
    u0 = apply_calcium_initial_conditions(settled_initial_state(params, gain, mode), calcium_ics)
    problem = ODEProblem(
        publication_driven_network_ode!,
        u0,
        (0.0, last(drive_time_s) * 1000.0),
        (params = sim_params, control_gain = gain, mode = mode, drive_t = drive_time_s, drive_v = drive_voltage_mv),
    )
    solution = solve(problem, RK4(); saveat = saveat_ms, reltol = 1e-6, abstol = 1e-6, maxiters = 10^9)
    states = Array(solution)
    solver_time_s = solution.t ./ 1000.0
    display_time_s = solver_time_s .* display_time_scale
    return DataFrame(
        mode = fill(mode_label(mode), length(display_time_s)),
        control_gain = fill(gain, length(display_time_s)),
        time_s = display_time_s,
        V0 = [interpolate_drive(t, drive_time_s, drive_voltage_mv) for t in solver_time_s],
        V1 = states[2, :],
        V3 = states[4, :],
        V4 = states[5, :],
    )
end

function simulate_scan_state_trace(params::ModelParams, gain::Float64, mode::ControlMode, config::SweepConfig; display_time_scale::Float64 = 1.0, calcium_ics = nothing)
    u0 = apply_calcium_initial_conditions(initial_state(), calcium_ics)
    sim = run_simulation(params, gain, mode, config; initial_u0 = u0)
    return DataFrame(
        mode = fill(mode_label(mode), length(sim.time_s)),
        control_gain = fill(gain, length(sim.time_s)),
        time_s = Float64.(sim.time_s) .* display_time_scale,
        V0 = Float64.(sim.V0),
        V1 = Float64.(sim.V1),
        V3 = Float64.(sim.V3),
        V4 = Float64.(sim.V4),
    )
end

function finite_cols(df::DataFrame, xcol::Symbol, ycol::Symbol)
    isempty(df) && return Float64[], Float64[]
    x = Float64.(df[!, xcol])
    y = Float64.(df[!, ycol])
    keep = isfinite.(x) .& isfinite.(y)
    return x[keep], y[keep]
end

function summary_vectors(summary::DataFrame, mode_name::String, ycol::Symbol)
    sub = sort(summary[summary.mode .== mode_name, :], :control_gain)
    return finite_cols(sub, :control_gain, ycol)
end

function raw_vectors(raw::DataFrame, mode_name::String, phase::String)
    sub = raw[(raw.mode .== mode_name) .& (raw.phase .== phase), :]
    return finite_cols(sub, :control_gain, :frequency_hz)
end
