using CSV
using DataFrames
using OrdinaryDiffEq

include(joinpath(@__DIR__, "legacy", "publication_plasticity_scan.jl"))

const CALIBRATION_OUTPUT_DIR = joinpath(@__DIR__, "outputs")
const CALIBRATION_PROVENANCE_CSV = joinpath(CALIBRATION_OUTPUT_DIR, "neuromod_calibration_provenance.csv")

const CITED_EXCITATORY_ALPHA = default_params().alphax
const PRESYNAPTIC_BASE_ALPHA = 0.011
const PRESYNAPTIC_ANCHOR_GAIN = CITED_EXCITATORY_ALPHA / PRESYNAPTIC_BASE_ALPHA - 1.0

const CITED_EXCITATORY_G = 0.002
const POSTSYNAPTIC_ANCHOR_GAIN = 3.0
const SATURATED_MOD_LEVEL = 1.0
const COMPENSATED_DIRECT_POST_BASE_G = CITED_EXCITATORY_G / (1.0 + POSTSYNAPTIC_ANCHOR_GAIN * SATURATED_MOD_LEVEL)
const CALIBRATED_X_SHIFT = -6.35
const PRESYNAPTIC_SI1_EXCITATORY_G = 0.0005
const POSTSYNAPTIC_SI1_EXCITATORY_G = CITED_EXCITATORY_G

function calibrated_params(; si1_excitatory_g::Float64 = POSTSYNAPTIC_SI1_EXCITATORY_G)
    updated_params(
        default_params();
        direct_post_base_g = COMPENSATED_DIRECT_POST_BASE_G,
        x_shift = CALIBRATED_X_SHIFT,
        x_shift1 = CALIBRATED_X_SHIFT,
        x_shift2 = CALIBRATED_X_SHIFT,
        x_shift3 = CALIBRATED_X_SHIFT,
        x_shift4 = CALIBRATED_X_SHIFT,
        g0 = si1_excitatory_g,
        t1_ms = 75_000.0,
        t2_ms = 1.0e12,
    )
end

function calibrated_params(mode::ControlMode)
    if mode == presynaptic
        return calibrated_params(; si1_excitatory_g = PRESYNAPTIC_SI1_EXCITATORY_G)
    elseif mode == postsynaptic
        return calibrated_params(; si1_excitatory_g = POSTSYNAPTIC_SI1_EXCITATORY_G)
    end
    error("Unsupported calibrated mode $(mode)")
end

function calibrated_scan_config()
    base = default_config()
    SweepConfig(
        900.0,
        base.saveat_ms,
        base.spike_threshold_mv,
        base.spike_refractory_s,
        base.burst_factor,
        base.onset_timeout_s,
        base.transient_bursts,
        base.measured_cycles,
    )
end

function calibrated_control_values(mode::ControlMode)
    if mode == presynaptic
        return sort!(unique(vcat(
            [0.0, PRESYNAPTIC_ANCHOR_GAIN],
            collect(range(0.30, 6.0; length = 72)),
        )))
    elseif mode == postsynaptic
        return sort!(unique(vcat(
            [0.0, 0.5, POSTSYNAPTIC_ANCHOR_GAIN],
            collect(range(0.75, 8.0; length = 64)),
        )))
    end
    error("Unsupported calibrated mode $(mode)")
end

function calibrated_spec(mode::ControlMode)
    label = mode == presynaptic ? "Presynaptic" : "Postsynaptic"
    ModeSweepSpec(mode, label, calibrated_control_values(mode), "Gain")
end

function append_calibration_columns!(summary::DataFrame, mode::ControlMode)
    mode_name = mode_label(mode)
    summary[!, :anchor_gain] .= mode == presynaptic ? PRESYNAPTIC_ANCHOR_GAIN : POSTSYNAPTIC_ANCHOR_GAIN
    summary[!, :anchor_target_alpha] .= mode == presynaptic ? CITED_EXCITATORY_ALPHA : NaN
    summary[!, :anchor_target_g] .= mode == postsynaptic ? CITED_EXCITATORY_G : NaN
    summary[!, :si1_excitatory_g] .= mode == presynaptic ? PRESYNAPTIC_SI1_EXCITATORY_G : POSTSYNAPTIC_SI1_EXCITATORY_G
    summary[!, :compensated_direct_post_base_g] .= mode == postsynaptic ? COMPENSATED_DIRECT_POST_BASE_G : NaN
    summary[!, :calibration_note] .= mode == postsynaptic ?
        "base g compensated so saturated anchor gain gives target effective g" :
        "anchor gain chosen so base alphax reaches cited alphax"
    summary[!, :mode] .= mode_name
    return summary
end

function calibration_provenance()
    DataFrame(
        key = [
            "cited_excitatory_alpha",
            "presynaptic_base_alpha",
            "presynaptic_anchor_gain",
            "cited_excitatory_g",
            "postsynaptic_anchor_gain",
            "saturated_mod_level_assumption",
            "compensated_direct_post_base_g",
            "calibrated_x_shift",
            "presynaptic_si1_excitatory_g",
            "postsynaptic_si1_excitatory_g",
        ],
        value = [
            CITED_EXCITATORY_ALPHA,
            PRESYNAPTIC_BASE_ALPHA,
            PRESYNAPTIC_ANCHOR_GAIN,
            CITED_EXCITATORY_G,
            POSTSYNAPTIC_ANCHOR_GAIN,
            SATURATED_MOD_LEVEL,
            COMPENSATED_DIRECT_POST_BASE_G,
            CALIBRATED_X_SHIFT,
            PRESYNAPTIC_SI1_EXCITATORY_G,
            POSTSYNAPTIC_SI1_EXCITATORY_G,
        ],
        note = [
            "target excitatory synaptic alpha used elsewhere in this model family",
            "presynaptic-mode unmodulated excitatory alpha",
            "CITED_EXCITATORY_ALPHA / PRESYNAPTIC_BASE_ALPHA - 1",
            "target effective excitatory conductance for common trace after suppressing no-input bursting",
            "postsynaptic common-trace gain",
            "sm / scale_sm is approximately 1 under sustained strong drive",
            "CITED_EXCITATORY_G / (1 + POSTSYNAPTIC_ANCHOR_GAIN * SATURATED_MOD_LEVEL)",
            "nearby excitability shift that suppresses no-input Si2 bursting in preliminary trace probes",
            "Si1 excitatory drive conductance chosen from preliminary presynaptic V1 scan probes after base network is quiet",
            "Si1 excitatory drive conductance kept at the cited conductance; gain 0 remains quiet in postsynaptic scan probes",
        ],
    )
end

function write_calibration_provenance()
    mkpath(CALIBRATION_OUTPUT_DIR)
    CSV.write(CALIBRATION_PROVENANCE_CSV, calibration_provenance())
end

function settled_initial_state(params::ModelParams, gain::Float64, mode::ControlMode; settle_s::Float64 = 300.0)
    drive_t = [0.0, settle_s]
    drive_v = [-44.0, -44.0]
    sim_params = mode_params(params, mode)
    problem = ODEProblem(
        publication_driven_network_ode!,
        initial_state(),
        (0.0, settle_s * 1000.0),
        (params = sim_params, control_gain = gain, mode = mode, drive_t = drive_t, drive_v = drive_v),
    )
    solution = solve(problem, RK4(); saveat = settle_s * 1000.0, reltol = 1e-6, abstol = 1e-6, maxiters = 10^9)
    return Vector(Array(solution)[:, end])
end
