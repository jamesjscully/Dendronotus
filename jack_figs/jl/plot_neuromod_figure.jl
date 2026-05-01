using CSV
using CairoMakie
using DataFrames

include(joinpath(@__DIR__, "neuromod_dashboard.jl"))

const NEUROMOD_CLEAN_FIGURE_PNG = joinpath(DASHBOARD_OUTPUT_DIR, "neuromod_clean_figure.png")

function main()
    raw = load_dashboard_raw()
    summary = load_dashboard_summary()
    fig = build_clean_figure(summary, raw, load_representative_traces())
    mkpath(DASHBOARD_OUTPUT_DIR)
    CairoMakie.save(NEUROMOD_CLEAN_FIGURE_PNG, fig)
    println("Saved $(NEUROMOD_CLEAN_FIGURE_PNG)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
