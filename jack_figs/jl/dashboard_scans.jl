using CSV
using DataFrames

const OUTPUT_DIR = joinpath(@__DIR__, "outputs")
const PRESYNAPTIC_POINTS_CSV = joinpath(OUTPUT_DIR, "presynaptic_scan_points.csv")
const PRESYNAPTIC_SUMMARY_CSV = joinpath(OUTPUT_DIR, "presynaptic_scan_summary.csv")
const POSTSYNAPTIC_POINTS_CSV = joinpath(OUTPUT_DIR, "postsynaptic_scan_points.csv")
const POSTSYNAPTIC_SUMMARY_CSV = joinpath(OUTPUT_DIR, "postsynaptic_scan_summary.csv")
const DASHBOARD_POINTS_CSV = joinpath(OUTPUT_DIR, "dashboard_scan_points.csv")
const DASHBOARD_SUMMARY_CSV = joinpath(OUTPUT_DIR, "dashboard_scan_summary.csv")

function main()
    isfile(PRESYNAPTIC_POINTS_CSV) || error("Missing $(PRESYNAPTIC_POINTS_CSV). Run presynaptic_scan.jl first.")
    isfile(PRESYNAPTIC_SUMMARY_CSV) || error("Missing $(PRESYNAPTIC_SUMMARY_CSV). Run presynaptic_scan.jl first.")
    isfile(POSTSYNAPTIC_POINTS_CSV) || error("Missing $(POSTSYNAPTIC_POINTS_CSV). Run postsynaptic_scan.jl first.")
    isfile(POSTSYNAPTIC_SUMMARY_CSV) || error("Missing $(POSTSYNAPTIC_SUMMARY_CSV). Run postsynaptic_scan.jl first.")

    points = vcat(
        CSV.read(PRESYNAPTIC_POINTS_CSV, DataFrame),
        CSV.read(POSTSYNAPTIC_POINTS_CSV, DataFrame),
    )
    summary = vcat(
        CSV.read(PRESYNAPTIC_SUMMARY_CSV, DataFrame),
        CSV.read(POSTSYNAPTIC_SUMMARY_CSV, DataFrame),
    )

    mkpath(OUTPUT_DIR)
    CSV.write(DASHBOARD_POINTS_CSV, points)
    CSV.write(DASHBOARD_SUMMARY_CSV, summary)
    println("Saved $(DASHBOARD_POINTS_CSV)")
    println("Saved $(DASHBOARD_SUMMARY_CSV)")
end

main()
