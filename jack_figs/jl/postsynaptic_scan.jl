using CSV
using DataFrames
using Base.Threads

include(joinpath(@__DIR__, "legacy", "publication_plasticity_scan.jl"))

const ROOT_OUTPUT_DIR = joinpath(@__DIR__, "outputs")
const POSTSYNAPTIC_POINTS_CSV = joinpath(ROOT_OUTPUT_DIR, "postsynaptic_scan_points.csv")
const POSTSYNAPTIC_SUMMARY_CSV = joinpath(ROOT_OUTPUT_DIR, "postsynaptic_scan_summary.csv")

const POSTSYNAPTIC_PARAMS = updated_params(
    default_params();
    direct_post_base_g = 0.00064,
    t1_ms = 75_000.0,
    t2_ms = 1.0e12,
)

const POSTSYNAPTIC_SPEC = ModeSweepSpec(
    postsynaptic,
    "Postsynaptic",
    vcat([0.0], collect(range(0.06060606060606061, 6.0; length = 99))),
    "Gain",
)

function rerun_with_12_threads_if_needed()
    if nthreads() == 12 || get(ENV, "NEUROMOD_ALREADY_RELAUNCHED", "0") == "1"
        return
    end
    cmd = `$(Base.julia_cmd()) --project=jack_figs --threads=12 $(abspath(@__FILE__))`
    println("Relaunching postsynaptic scan with 12 threads...")
    run(addenv(cmd, "NEUROMOD_ALREADY_RELAUNCHED" => "1"))
    exit()
end

function generate_postsynaptic_data()
    gains = POSTSYNAPTIC_SPEC.control_values
    raw_parts = Vector{DataFrame}(undef, length(gains))
    summary_parts = Vector{DataFrame}(undef, length(gains))
    config = default_config()

    @threads for idx in eachindex(gains)
        gain = gains[idx]
        sim = run_simulation(POSTSYNAPTIC_PARAMS, gain, POSTSYNAPTIC_SPEC.mode, config)
        run_raw = publication_raw_points(sim, gain, POSTSYNAPTIC_SPEC.mode, POSTSYNAPTIC_PARAMS, config)
        raw_parts[idx] = run_raw
        summary_parts[idx] = publication_summary_row(sim, run_raw, gain, POSTSYNAPTIC_SPEC.mode)
    end

    return vcat(raw_parts...), vcat(summary_parts...)
end

function main()
    rerun_with_12_threads_if_needed()
    println("Running postsynaptic scan with $(nthreads()) threads...")
    raw, summary = generate_postsynaptic_data()
    mkpath(ROOT_OUTPUT_DIR)
    CSV.write(POSTSYNAPTIC_POINTS_CSV, raw)
    CSV.write(POSTSYNAPTIC_SUMMARY_CSV, summary)
    println("Saved $(POSTSYNAPTIC_POINTS_CSV)")
    println("Saved $(POSTSYNAPTIC_SUMMARY_CSV)")
end

main()
