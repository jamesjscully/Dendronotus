using CSV
using DataFrames
using Statistics

include(joinpath(@__DIR__, "compare_synaptic_plasticity_models.jl"))

const PUB_SCAN_RAW_CSV = joinpath(OUTPUT_DIR, "plasticity_publication_points.csv")
const PUB_SCAN_SUMMARY_CSV = joinpath(OUTPUT_DIR, "plasticity_publication_summary.csv")

function publication_params()
    updated_params(
        default_params();
        direct_post_base_g = 0.00064,
    )
end

function publication_specs()
    [
        ModeSweepSpec(
            presynaptic,
            "Presynaptic",
            vcat([0.0], collect(range(0.3050847457627119, 6.0; length = 119))),
            "Gain",
        ),
        ModeSweepSpec(
            postsynaptic,
            "Postsynaptic",
            vcat([0.0], collect(range(0.06060606060606061, 6.0; length = 99))),
            "Gain",
        ),
    ]
end

function empty_publication_raw()
    DataFrame(
        mode = String[],
        control_gain = Float64[],
        phase = String[],
        interval_start_s = Float64[],
        frequency_hz = Float64[],
        driver_frequency_hz = Float64[],
    )
end

function empty_publication_summary()
    DataFrame(
        mode = String[],
        control_gain = Float64[],
        onset_failed = Bool[],
        n_pre_intervals = Int[],
        mean_pre_frequency_hz = Float64[],
        mean_pre_driver_frequency_hz = Float64[],
        n_post_intervals = Int[],
        mean_post_frequency_hz = Float64[],
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
        if !isfinite(v_prev) || !isfinite(v_now) || !isfinite(t_prev) || !isfinite(t_now)
            continue
        end
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

function interval_spike_frequency(spike_times_s::AbstractVector{<:Real}, t_start::Float64, t_stop::Float64)
    duration = t_stop - t_start
    duration <= 0 && return NaN
    spike_count = Base.count(t -> t >= t_start && t < t_stop, spike_times_s)
    return spike_count / duration
end

function pre_interval_metrics(sim, control_gain::Float64, mode::ControlMode, config::SweepConfig)
    raw = empty_publication_raw()
    required_bursts = config.transient_bursts + config.measured_cycles
    bursts = detect_bursts_from_spikes(sim.spike_times_s, config)
    complete_start_times = complete_burst_start_times(bursts)
    if sim.onset_failed || length(complete_start_times) < required_bursts
        return raw
    end

    measured_start_times = complete_start_times[(config.transient_bursts + 1):(config.transient_bursts + config.measured_cycles)]
    driver_spikes = threshold_crossing_spikes(sim.time_s, sim.V0)

    for k in 1:(length(measured_start_times) - 1)
        interval_start = measured_start_times[k]
        interval_stop = measured_start_times[k + 1]
        interval_stop <= interval_start && continue
        push!(raw, (
            mode_label(mode),
            control_gain,
            "pre",
            interval_start,
            1.0 / (interval_stop - interval_start),
            interval_spike_frequency(driver_spikes, interval_start, interval_stop),
        ))
    end

    return raw
end

function post_withdrawal_intervals(sim, config::SweepConfig)
    if sim.onset_failed || isnan(sim.stimulus_off_s) || isnan(sim.withdrawal_end_s)
        return Float64[], Float64[]
    end

    bursts = detect_bursts_from_spikes(sim.spike_times_s, config)
    starts = complete_burst_start_times(bursts)
    ends = complete_burst_end_times(bursts)
    keep = findall(i -> starts[i] >= sim.stimulus_off_s && ends[i] <= sim.withdrawal_end_s, eachindex(starts))
    starts = starts[keep]
    if length(starts) < 2
        return Float64[], Float64[]
    end

    interval_starts = starts[1:(end - 1)]
    freqs = 1.0 ./ diff(starts)
    return interval_starts, freqs
end

function publication_raw_points(sim, control_gain::Float64, mode::ControlMode, params::ModelParams, config::SweepConfig)
    raw = empty_publication_raw()

    append!(raw, pre_interval_metrics(sim, control_gain, mode, config))

    post_starts, post_freqs = post_withdrawal_intervals(sim, config)
    for (t, f) in zip(post_starts, post_freqs)
        push!(raw, (
            mode_label(mode),
            control_gain,
            "post",
            t,
            f,
            NaN,
        ))
    end

    return raw
end

function publication_summary_row(sim, raw::DataFrame, control_gain::Float64, mode::ControlMode)
    pre = raw[raw.phase .== "pre", :]
    post = raw[raw.phase .== "post", :]

    mean_pre = isempty(pre) ? 0.0 : mean(pre.frequency_hz)
    mean_pre_driver = isempty(pre) ? 0.0 : mean(skipmissing(pre.driver_frequency_hz))
    mean_post = isempty(post) ? 0.0 : mean(post.frequency_hz)

    return DataFrame(
        mode = [mode_label(mode)],
        control_gain = [control_gain],
        onset_failed = [sim.onset_failed],
        n_pre_intervals = [nrow(pre)],
        mean_pre_frequency_hz = [mean_pre],
        mean_pre_driver_frequency_hz = [mean_pre_driver],
        n_post_intervals = [nrow(post)],
        mean_post_frequency_hz = [mean_post],
    )
end

function generate_publication_data(config::SweepConfig = default_config(), params::ModelParams = publication_params(), specs::Vector{ModeSweepSpec} = publication_specs())
    raw = empty_publication_raw()
    summary = empty_publication_summary()

    for spec in specs
        println("Running publication $(spec.label) sweep...")
        for control_gain in spec.control_values
            println("  control_gain = $(control_gain)")
            sim = run_simulation(params, control_gain, spec.mode, config)
            run_raw = publication_raw_points(sim, control_gain, spec.mode, params, config)
            append!(raw, run_raw)
            append!(summary, publication_summary_row(sim, run_raw, control_gain, spec.mode))
        end
    end

    return raw, summary
end

function save_publication_data(raw::DataFrame, summary::DataFrame)
    mkpath(OUTPUT_DIR)
    CSV.write(PUB_SCAN_RAW_CSV, raw)
    CSV.write(PUB_SCAN_SUMMARY_CSV, summary)
end

function load_publication_raw()
    CSV.read(PUB_SCAN_RAW_CSV, DataFrame)
end

function load_publication_summary()
    CSV.read(PUB_SCAN_SUMMARY_CSV, DataFrame)
end

function main(args = ARGS)
    mode = isempty(args) ? "both" : lowercase(args[1])

    if mode == "generate" || mode == "both"
        raw, summary = generate_publication_data()
        save_publication_data(raw, summary)
        println("Saved:")
        println("  $(PUB_SCAN_RAW_CSV)")
        println("  $(PUB_SCAN_SUMMARY_CSV)")
        return
    end

    if mode == "inspect"
        summary = load_publication_summary()
        show(first(summary, 20), allrows = true, allcols = true)
        println()
        return
    end

    error("Unknown mode `$(mode)`. Use `generate`, `both`, or `inspect`.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
