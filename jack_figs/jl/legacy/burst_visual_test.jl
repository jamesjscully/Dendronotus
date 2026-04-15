using Plots
include("neural_simulation.jl")

# Test parameters - you can modify these to test different conditions
test_gm_values = [0.0, 0.01, 0.05]  # Test a few different gm values
tmax = 200000.0
tsamp = 2.0

# Burst detection parameters
burst_threshold = -20.0
min_interburst_interval = 5000.0
analysis_start_time = 50.0
analysis_end_time = 190.0

println("=== BURST DETECTION VISUAL TEST ===")
println("Creating visual validation plots for burst detection algorithm")
println("Please manually count the bursts and compare with algorithm results")

for (idx, gm) in enumerate(test_gm_values)
    println("\n--- Testing gm = $gm ---")
    
    # Run simulation
    time, vv0, vv1, vv2, vv3, vv4, sol = run_neural_simulation(gm, tmax=tmax, tsamp=tsamp)
    
    # Detect bursts
    burst_frequency, num_bursts, peak_locs, peaks = detect_bursts(
        time, vv4, 
        burst_threshold=burst_threshold,
        min_interburst_interval=min_interburst_interval,
        analysis_start_time=analysis_start_time,
        analysis_end_time=analysis_end_time,
        tsamp=tsamp
    )
    
    # Create comprehensive test plot
    p_test = plot(time, vv4, 
                 linewidth=1.5, 
                 linecolor=:blue,
                 xlabel="Time (s)",
                 ylabel="Voltage (mV)",
                 title="Burst Detection Test: gm = $gm\nAlgorithm Count = $num_bursts bursts, Frequency = $(round(burst_frequency, digits=3)) Hz",
                 legend=:topright,
                 size=(1000, 600),
                 label="V4 voltage trace")
    
    # Add horizontal line at burst threshold
    hline!(p_test, [burst_threshold], 
           linewidth=2, 
           linestyle=:dash, 
           linecolor=:red,
           label="Burst threshold (-20 mV)")
    
    # Overlay detected peaks
    if length(peak_locs) > 0
        scatter!(p_test, time[peak_locs], vv4[peak_locs],
                markersize=8,
                markercolor=:red,
                markerstrokewidth=2,
                markerstrokecolor=:darkred,
                label="Detected peaks ($(length(peak_locs)) total)")
    end
    
    # Add analysis window markers
    vline!(p_test, [analysis_start_time, analysis_end_time], 
           linewidth=3, 
           linestyle=:dot, 
           linecolor=:green,
           alpha=0.8,
           label="Analysis window")
    
    # Add shaded analysis region
    analysis_y_min, analysis_y_max = ylims(p_test)
    plot!(p_test, [analysis_start_time, analysis_end_time, analysis_end_time, analysis_start_time, analysis_start_time], 
          [analysis_y_min, analysis_y_min, analysis_y_max, analysis_y_max, analysis_y_min],
          fillalpha=0.1, 
          fillcolor=:green,
          linewidth=0,
          label="")
    
    # Count peaks within analysis window for verification
    analysis_peak_indices = findall(x -> x >= analysis_start_time && x <= analysis_end_time, time[peak_locs])
    peaks_in_window = length(analysis_peak_indices)
    
    # Add text annotations
    annotate!(p_test, [(analysis_start_time + 10, maximum(vv4) * 0.9, 
                       text("Analysis Window\n$(analysis_start_time)s - $(analysis_end_time)s\n$peaks_in_window bursts counted", 
                            :left, 10, :darkgreen))])
    
    # Highlight peaks within analysis window
    if length(analysis_peak_indices) > 0
        analysis_peak_locs = peak_locs[analysis_peak_indices]
        scatter!(p_test, time[analysis_peak_locs], vv4[analysis_peak_locs],
                markersize=10,
                markercolor=:orange,
                markerstrokewidth=3,
                markerstrokecolor=:darkorange,
                markershape=:star5,
                label="Counted bursts ($peaks_in_window)")
    end
    
    display(p_test)
    
    # Save the test plot
    filename = "burst_visual_test_gm_$(gm).png"
    savefig(p_test, filename)
    println("Visual test plot saved as: $filename")
    
    # Print detailed analysis
    println("DETAILED ANALYSIS:")
    println("  Total peaks detected: $(length(peak_locs))")
    println("  Peaks in analysis window ($analysis_start_time-$(analysis_end_time)s): $peaks_in_window")
    println("  Algorithm burst count: $num_bursts")
    println("  Calculated frequency: $(round(burst_frequency, digits=3)) Hz")
    println("  Analysis duration: $(analysis_end_time - analysis_start_time) seconds")
    
    if length(peak_locs) > 0
        println("  Peak times (all): ", round.(time[peak_locs], digits=1))
        if length(analysis_peak_indices) > 0
            analysis_peak_times = time[peak_locs[analysis_peak_indices]]
            println("  Peak times (in window): ", round.(analysis_peak_times, digits=1))
        end
    end
    
    println("\n>>> MANUAL VERIFICATION INSTRUCTIONS:")
    println("    1. Look at the plot and manually count voltage peaks above the red dashed line")
    println("    2. Only count peaks between the green dotted vertical lines (analysis window)")
    println("    3. Orange stars show the peaks the algorithm counted")
    println("    4. Compare your manual count with the algorithm count: $num_bursts")
    println("    5. Check that orange stars align with actual voltage peaks")
end

println("\n=== VISUAL TEST COMPLETE ===")
println("Please review the generated plots and verify the burst detection accuracy.")
println("The algorithm should correctly identify voltage peaks above threshold within the analysis window.")