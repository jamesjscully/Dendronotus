# Workflow:
#   From the repo root, launch the dashboard with:
#       julia --project=jack_figs jack_figs/jl/neuromod_dashboard.jl
#
#   For interactive development, keep Julia open and reload only the layer you are editing:
#       julia --project=jack_figs -i
#       include("jack_figs/jl/neuromod_dashboard_numerics.jl")
#       include("jack_figs/jl/neuromod_dashboard_plotting.jl")
#       dashboard = build_dashboard()
#       screen = display(dashboard.fig)
#
#   When changing only widgets/layout/plot behavior, rerun:
#       include("jack_figs/jl/neuromod_dashboard_plotting.jl")
#       dashboard = build_dashboard()
#       screen = display(dashboard.fig)
#
#   The live handle exposes observables and actions, e.g.:
#       GLMakie.set_close_to!(dashboard.param_controls.sliders[1], -6.25)
#       # trace slider: gain
#       # parameter sliders: x_shift1/2, x_shift3/4, g0 presyn, g0 postsyn,
#       # alpha1, beta1, alpham, betam, Ca_shift1/2, Ca_shift3/4,
#       # g41/g32 base, alphax, t1
#       dashboard.selected[] = ("postsynaptic", 3.0)
#       # Trace panels update reactively from slider changes; call this only
#       # when changing non-slider state directly from the REPL.
#       dashboard.actions.request_trace_update!()
#       dashboard.actions.run_scan!()
#
#   Use `julia --project=jack_figs jack_figs/jl/neuromod_dashboard.jl clean`
#   or the dashboard save button to write a PNG without UI widgets.

include(joinpath(@__DIR__, "neuromod_dashboard_numerics.jl"))
include(joinpath(@__DIR__, "neuromod_dashboard_plotting.jl"))

function main(args = ARGS)
    mode = isempty(args) ? "interactive" : lowercase(args[1])
    if mode == "interactive"
        dashboard = build_dashboard()
        screen = display(dashboard.fig)
        wait(screen)
        return dashboard
    elseif mode == "clean"
        clean = build_clean_figure(load_dashboard_summary(), load_dashboard_raw(), load_representative_traces())
        mkpath(DASHBOARD_OUTPUT_DIR)
        CairoMakie.save(DASHBOARD_CLEAN_PNG, clean)
        println("Saved $(DASHBOARD_CLEAN_PNG)")
        return clean
    end
    error("Unknown mode `$(mode)`. Use `interactive` or `clean`.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
