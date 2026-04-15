using Plots, CSV, DataFrames, Statistics
include("neural_simulation.jl")

# Parameter sweep setup
gm_values = [0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1] # Different gm levels to test
burst_frequencies = zeros(length(gm_values)) # Store burst frequencies

# Simulation parameters
tmax = 200.0 # max time in seconds
tsamp = 2.0  # sampling interval to save data in seconds

# Burst detection parameters
burst_threshold = -20.0 # Voltage threshold for detecting spikes/bursts
min_burst_duration = 100.0 # Minimum duration in msec to count as burst
min_interburst_interval = 500.0 # Minimum time between bursts in msec

# Loop over different gm values
println("Running simulations for different gm values...")
for (gm_idx, gm) in enumerate(gm_values)
    println("Running simulation $gm_idx/$(length(gm_values)): gm = $gm")
    
    # Run simulation
    time, vv0, vv1, vv2, vv3, vv4, sol = run_neural_simulation(gm, tmax=tmax, tsamp=tsamp)
    
    # Detect bursts using cell V4 as representative
    analysis_start_time = 50.0 # Start analysis after 50 seconds to avoid initial transients
    analysis_end_time = tmax - 10.0 # End 10 seconds before simulation end
    
    burst_frequency, num_bursts, peak_locs, peaks = detect_bursts(
        time, vv4, 
        burst_threshold=burst_threshold,
        min_interburst_interval=min_interburst_interval,
        analysis_start_time=analysis_start_time,
        analysis_end_time=analysis_end_time,
        tsamp=tsamp
    )
    
    burst_frequencies[gm_idx] = burst_frequency
    
    println("  Detected $num_bursts bursts, frequency = $(round(burst_frequency, digits=3)) Hz")
end

# Create scatterplot
p = scatter(gm_values, burst_frequencies, 
           markersize=8, 
           markerstrokewidth=2,
           markercolor=:steelblue,
           markerstrokecolor=:black,
           xlabel="gm Parameter Value",
           ylabel="Burst Frequency (Hz)",
           title="Effect of gm Parameter on Network Burst Frequency",
           legend=false,
           grid=true,
           gridwidth=1,
           gridcolor=:lightgray,
           framestyle=:box,
           size=(800, 600))

# Add trend line if there's a clear relationship
if length(gm_values) > 2
    # Linear fit
    A = [ones(length(gm_values)) gm_values]
    coeffs = A \ burst_frequencies
    trend_line = coeffs[1] .+ coeffs[2] .* gm_values
    plot!(p, gm_values, trend_line, 
          linewidth=2, 
          linestyle=:dash, 
          linecolor=:red)
    
    # Calculate correlation
    correlation = cor(gm_values, burst_frequencies)
    annotate!(p, [(0.05*maximum(gm_values), 0.95*maximum(burst_frequencies), 
                   text("Correlation: $(round(correlation, digits=3))", 
                        :left, 10, :black))])
end

display(p)

# Display results table
println("\n=== RESULTS SUMMARY ===")
println("gm Value\tBurst Frequency (Hz)")
println("--------\t----------------")
for i in 1:length(gm_values)
    println("$(round(gm_values[i], digits=4))\t\t$(round(burst_frequencies[i], digits=4))")
end

# Save results
results_df = DataFrame(gm_parameter=gm_values, burst_frequency_Hz=burst_frequencies)
CSV.write("gm_burst_analysis_results.csv", results_df)
println("\nResults saved to: gm_burst_analysis_results.csv")

# Save the main plot
savefig(p, "gm_parameter_sweep.png")
println("Parameter sweep plot saved as: gm_parameter_sweep.png")