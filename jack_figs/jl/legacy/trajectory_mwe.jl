using Plots
using Plots.PlotMeasures

include(joinpath(@__DIR__, "compare_synaptic_plasticity_models.jl"))

function mwe_config()
    base = default_config()
    SweepConfig(
        120.0,
        250.0,
        base.spike_threshold_mv,
        base.spike_refractory_s,
        base.burst_factor,
        base.transient_bursts,
        base.measured_cycles,
    )
end

function main(args = ARGS)
    length(args) >= 2 || error("Usage: trajectory_mwe.jl <mode> <control_gain>")
    mode = parse_mode_label(args[1])
    control_gain = parse(Float64, args[2])

    config = mwe_config()
    params = default_params()
    sim = run_simulation(params, control_gain, mode, config)
    bursts = detect_bursts_from_spikes(sim.spike_times_s, config)
    traces = effective_excitation_trace(sim.sm, sim.s4, control_gain, mode, sim.params)

    complete_starts = complete_burst_start_times(bursts)
    incomplete_starts = incomplete_burst_start_times(bursts)

    p1 = plot(
        sim.time_s,
        sim.V1;
        label = "V1",
        color = :forestgreen,
        linewidth = 1.2,
        xlabel = "Time (s)",
        ylabel = "Voltage (mV)",
        title = "Trajectory MWE: $(mode_label(mode)) gain=$(control_gain)",
        framestyle = :box,
        legend = :topright,
    )
    plot!(p1, sim.time_s, sim.V4; label = "V4", color = :firebrick, linewidth = 1.2)
    for t in sim.spike_times_s
        vline!(p1, [t]; color = RGBA(0.5, 0.5, 0.5, 0.18), linewidth = 0.8, label = false)
    end
    for t in complete_starts
        vline!(p1, [t]; color = :black, linewidth = 1.4, label = false)
    end
    for t in incomplete_starts
        vline!(p1, [t]; color = :orange, linewidth = 1.8, label = false)
    end

    p2 = plot(
        sim.time_s,
        traces.strength;
        label = "g_eff * s4",
        color = :purple,
        linewidth = 1.4,
        xlabel = "Time (s)",
        ylabel = "Excitation",
        framestyle = :box,
        legend = :topright,
    )
    plot!(p2, sim.time_s, sim.sm; label = "sm", color = :darkorange, linewidth = 1.2)
    plot!(p2, sim.time_s, traces.gain; label = "gain", color = :teal, linewidth = 1.2)

    fig = plot(
        p1,
        p2;
        layout = (2, 1),
        size = (1800, 900),
        left_margin = 6mm,
        right_margin = 4mm,
        bottom_margin = 8mm,
        top_margin = 6mm,
    )

    outfile = joinpath(OUTPUT_DIR, "trajectory_mwe_$(mode_label(mode))_$(replace(string(control_gain), '.' => 'p')).png")
    mkpath(dirname(outfile))
    savefig(fig, outfile)

    println("mode=$(mode_label(mode)) gain=$(control_gain)")
    println("spikes=$(length(sim.spike_times_s))")
    println("complete_bursts=$(count(identity, bursts.is_complete))")
    println("incomplete_bursts=$(count(!, bursts.is_complete))")
    println("burst_counts=$(bursts.spike_counts)")
    println("saved=$(outfile)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
