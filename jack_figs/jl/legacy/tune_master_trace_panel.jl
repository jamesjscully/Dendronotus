ENV["NEUROMOD_MASTER_NOAUTORUN"] = "1"
include(joinpath(@__DIR__, "neuromod_master_figure.jl"))

function burst_metrics(time_s::AbstractVector{<:Real}, voltage_mv::AbstractVector{<:Real})
    starts, _ = detect_bursts_from_voltage(time_s, voltage_mv; threshold_mv = -40.0, min_down_time_s = 1.0)
    burst_times = Float64.(time_s[starts])
    if length(burst_times) < 2
        return (n = length(burst_times), onset = Inf, mean_freq = 0.0)
    end
    freqs = 1.0 ./ diff(burst_times)
    return (n = length(burst_times), onset = burst_times[1], mean_freq = mean(freqs))
end

function score_candidate(metrics, target_onset::Float64, target_freq::Float64)
    onset_term = metrics.onset == Inf ? 1000.0 : abs(metrics.onset - target_onset)
    freq_term = abs(metrics.mean_freq - target_freq)
    sparse_penalty = metrics.n < 3 ? 10.0 : 0.0
    return onset_term + 40.0 * freq_term + sparse_penalty
end

function run_tuning()
    bio_t1, bio_v1 = load_two_column_txt(BIO_SI1_TXT)
    bio_t2, bio_v2 = load_two_column_txt(BIO_SI2_TXT)
    bio_t1 = bio_t1 .- first(bio_t1)
    bio_t2 = bio_t2 .- first(bio_t2)
    target = burst_metrics(bio_t2, bio_v2)
    println("Biological target: onset=$(round(target.onset, digits=3)) s, mean_freq=$(round(target.mean_freq, digits=4)) Hz, n=$(target.n)")

    pre_candidates = NamedTuple[]
    for alpha_scale in (1.0, 1.25, 1.5, 1.75, 2.0)
        for gm_scale in (1.0, 1.25, 1.5, 1.75)
            for g_scale in (1.0, 1.25, 1.5)
                trace = simulate_compare_model_driven("presynaptic", bio_t1, bio_v1; alpha_scale = alpha_scale, gm_scale = gm_scale, g_scale = g_scale, step_ms = 0.2, saveat_ms = 5.0)
                metrics = burst_metrics(trace[!, :time_s], trace[!, :V1])
                push!(pre_candidates, (alpha_scale = alpha_scale, gm_scale = gm_scale, g_scale = g_scale, onset = metrics.onset, mean_freq = metrics.mean_freq, n = metrics.n, score = score_candidate(metrics, target.onset, target.mean_freq)))
            end
        end
    end

    post_candidates = NamedTuple[]
    for alpha_scale in (1.0, 1.25, 1.5)
        for gm_scale in (1.0, 1.25, 1.5, 1.75, 2.0)
            for g_scale in (1.0, 1.25, 1.5, 2.0, 2.5)
                trace = simulate_compare_model_driven("postsynaptic", bio_t1, bio_v1; alpha_scale = alpha_scale, gm_scale = gm_scale, g_scale = g_scale, step_ms = 0.2, saveat_ms = 5.0)
                metrics = burst_metrics(trace[!, :time_s], trace[!, :V1])
                push!(post_candidates, (alpha_scale = alpha_scale, gm_scale = gm_scale, g_scale = g_scale, onset = metrics.onset, mean_freq = metrics.mean_freq, n = metrics.n, score = score_candidate(metrics, target.onset, target.mean_freq)))
            end
        end
    end

    sort!(pre_candidates, by = x -> x.score)
    sort!(post_candidates, by = x -> x.score)

    println("\nBest presynaptic candidates:")
    for row in first(pre_candidates, min(10, length(pre_candidates)))
        println(row)
    end

    println("\nBest postsynaptic candidates:")
    for row in first(post_candidates, min(10, length(post_candidates)))
        println(row)
    end
end

run_tuning()
