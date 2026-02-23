include("LaNMM_Engine.jl")
include("LaNMM_Analyzer.jl")

using Dates
using Serialization
using Printf

# ==============================================================================
# 1. USER CONFIGURATION
# ==============================================================================
const intrinsic_params = Dict(
    "condition" => "healthy",
    "v0_default" => 6.0,
    "v0_p2" => 1.0,
    "fmax" => 5.0,
    "r_slope" => 0.56
)

const job_params = Dict(
    "job_title" => "julia_lanmm_sweep",
    "mu_p1_values" => collect(50.0:10.0:400.0),
    "mu_p2_values" => collect(50.0:10.0:400.0),
    "tmax" => 100.0,   # change to 300.0 for production sweeps
    "dt" => 0.001,
    "discard" => 1.0,
    "alpha_band" => (8.0, 12.0),
    "gamma_band" => (30.0, 50.0),
    "quiet_progress" => false # true => fewer progress updates
)

# ==============================================================================
# 2. JOB RUNNER
# ==============================================================================
"""
Run a parameter sweep and generate the same plot families as the Python pipeline:
- single-run diagnostics (`inputs`, `v`, `u_external`, `psd`)
- coupling heatmaps (`SEC`, `EEC`)
- band-power heatmaps (`P1/P2` x `alpha/gamma`)
- PEIX heatmaps (`P1`, `P2`)
- peak-frequency heatmaps (`alpha`, `gamma`)
"""
function run_sweep_job(intrinsic, job)
    job_title = job["job_title"]
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    output_dir = "$(job_title)_$(timestamp)"
    mkpath(output_dir)

    println("Saving parameters to $output_dir/parameters.txt")
    open(joinpath(output_dir, "parameters.txt"), "w") do f
        println(f, "# Julia LaNMM sweep parameters")
        println(f, "intrinsic = ", intrinsic)
        println(f, "job = ", job)
    end

    # One nominal run for detailed per-simulation plots.
    nominal = run_unified_simulation(
        condition=intrinsic["condition"],
        v0_default=intrinsic["v0_default"],
        v0_p2=intrinsic["v0_p2"],
        fmax=intrinsic["fmax"],
        r_slope=intrinsic["r_slope"],
        mu_e1=270.0,
        mu_e2=90.0,
        tmax=job["tmax"],
        dt=job["dt"],
        discard=job["discard"]
    )
    plot_sim_results(nominal; save_dir=output_dir, xlims=(5, 10))

    m1_vals = job["mu_p1_values"]
    m2_vals = job["mu_p2_values"]

    println("Running sweep for $(length(m1_vals) * length(m2_vals)) parameter pairs...")
    sweep_results = Dict{Tuple{Float64,Float64},Any}()
    total_pairs = length(m1_vals) * length(m2_vals)
    done_pairs = 0
    t0 = time()
    quiet_progress = get(job, "quiet_progress", false)
    # Verbose mode: ~0.5% increments. Quiet mode: ~10% increments.
    update_every = quiet_progress ? max(1, total_pairs ÷ 10) : max(1, total_pairs ÷ 200)

    function fmt_seconds(sec::Real)
        sec_i = max(0, round(Int, sec))
        h = sec_i ÷ 3600
        m = (sec_i % 3600) ÷ 60
        s = sec_i % 60
        if h > 0
            return @sprintf("%dh%02dm%02ds", h, m, s)
        elseif m > 0
            return @sprintf("%dm%02ds", m, s)
        else
            return @sprintf("%ds", s)
        end
    end

    function fmt_finish_time(unix_ts::Real)
        dt_finish = Dates.unix2datetime(floor(Int, unix_ts))
        return Dates.format(dt_finish, "HH:MM:SS")
    end

    function print_progress(done::Int, total::Int, t_start::Float64; quiet::Bool=false)
        frac = done / total
        pct = 100 * frac
        elapsed = time() - t_start
        rate = done / max(elapsed, 1e-9)
        eta = (total - done) / max(rate, 1e-9)
        avg_per_sim = elapsed / max(done, 1)
        finish_clock = fmt_finish_time(time() + eta)

        width = 30
        filled = clamp(round(Int, width * frac), 0, width)
        bar = repeat("=", filled) * repeat(".", width - filled)

        if quiet
            @printf(
                "Sweep progress [%s] %6.2f%% (%d/%d) | elapsed %s | ETA %s | avg %.3fs/sim | finish ~%s\n",
                bar, pct, done, total, fmt_seconds(elapsed), fmt_seconds(eta), avg_per_sim, finish_clock
            )
        else
            @printf(
                "\rSweep progress [%s] %6.2f%% (%d/%d) | elapsed %s | ETA %s | avg %.3fs/sim | finish ~%s",
                bar, pct, done, total, fmt_seconds(elapsed), fmt_seconds(eta), avg_per_sim, finish_clock
            )
            flush(stdout)
            if done == total
                println()
            end
        end
    end

    for m1 in m1_vals, m2 in m2_vals
        sweep_results[(m1, m2)] = run_unified_simulation(
            condition=intrinsic["condition"],
            v0_default=intrinsic["v0_default"],
            v0_p2=intrinsic["v0_p2"],
            fmax=intrinsic["fmax"],
            r_slope=intrinsic["r_slope"],
            mu_e1=m1,
            mu_e2=m2,
            tmax=job["tmax"],
            dt=job["dt"],
            discard=job["discard"]
        )
        done_pairs += 1
        if done_pairs == 1 || done_pairs % update_every == 0 || done_pairs == total_pairs
            print_progress(done_pairs, total_pairs, t0; quiet=quiet_progress)
        end
    end

    println("Computing analysis products...")
    couplings = analyze_sweep_couplings(
        sweep_results;
        alpha_band=job["alpha_band"],
        gamma_band=job["gamma_band"]
    )
    power_results = analyze_sweep_power(
        sweep_results;
        alpha_band=job["alpha_band"],
        gamma_band=job["gamma_band"],
        method=:bandpass_hilbert
    )
    peix_results = sweep_peix(
        sweep_results;
        r_slope=intrinsic["r_slope"],
        v0_default=intrinsic["v0_default"],
        v0_p2=intrinsic["v0_p2"]
    )
    freq_data = analyze_peak_frequencies(sweep_results; fs=1.0 / job["dt"])

    println("Saving plots...")
    plot_coupling_heatmaps(
        couplings;
        title_str="Couplings: $(job_title)",
        save_path=joinpath(output_dir, "couplings.png")
    )
    plot_power_heatmaps_bottom_cbar(
        power_results;
        title_str="Band Power: $(job_title)",
        save_path=joinpath(output_dir, "power.png"),
        log_scale=true
    )
    plot_peix_heatmaps(
        peix_results;
        title_str="PEIX: $(job_title)",
        save_path=joinpath(output_dir, "peix.png")
    )
    plot_frequency_heatmaps(
        freq_data;
        save_path=joinpath(output_dir, "frequency_heatmaps.png")
    )

    println("Saving raw arrays...")
    serialize(joinpath(output_dir, "sweep_results.jls"), sweep_results)
    serialize(
        joinpath(output_dir, "analysis_results.jls"),
        (couplings=couplings, power_results=power_results, peix_results=peix_results, freq_data=freq_data)
    )

    println("Job completed successfully. Outputs in: $output_dir")
    return (output_dir=output_dir, sweep_results=sweep_results)
end

# ==============================================================================
# 3. EXECUTE
# ==============================================================================
@time run_sweep_job(intrinsic_params, job_params)