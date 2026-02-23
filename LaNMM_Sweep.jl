using Dates
using Serialization
using Printf
using TOML
if !isdefined(@__MODULE__, :IntrinsicConfig)
    include("LaNMM_Types.jl")
end
if !isdefined(@__MODULE__, :get_column_parameters)
    include("LaNMM_parameters.jl")
end
if !isdefined(@__MODULE__, :run_unified_simulation)
    include("LaNMM_Engine.jl")
end
if !isdefined(@__MODULE__, :analyze_sweep_couplings)
    include("LaNMM_Analyzer.jl")
end
if !isdefined(@__MODULE__, :plot_sim_results)
    include("LaNMM_Plots.jl")
end

"""
Public API (this file)
----------------------
- `collect_run_parameters(intrinsic, job)`
- `save_parameters_snapshot(path, intrinsic, job)`
- `validate_configs(intrinsic, job)`
- `run_sweep_job(intrinsic, job)`

This file orchestrates a complete sweep workflow and should not define ODE or
analysis equations directly.
"""

"""
Pretty-print a duration in seconds as `HhMMmSSs`, `MmSSs`, or `Ss`.
"""
function _fmt_seconds(sec::Real)
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

"""
Format unix timestamp to clock time string.
"""
function _fmt_finish_time(unix_ts::Real)
    dt_finish = Dates.unix2datetime(floor(Int, unix_ts))
    return Dates.format(dt_finish, "HH:MM:SS")
end

function _to_toml_value(x)
    if x isa Dict
        out = Dict{String,Any}()
        for (k, v) in x
            out[string(k)] = _to_toml_value(v)
        end
        return out
    elseif x isa NamedTuple
        out = Dict{String,Any}()
        for name in keys(x)
            out[string(name)] = _to_toml_value(getfield(x, name))
        end
        return out
    elseif x isa Tuple
        return [_to_toml_value(v) for v in x]
    elseif x isa AbstractVector
        return [_to_toml_value(v) for v in x]
    elseif x isa Symbol
        return string(x)
    else
        return x
    end
end

"""
Build a full parameter snapshot for reproducibility.

Includes:
- user-provided intrinsic/sweep configs
- effective column-level parameters
- derived synapse arrays used by ODE
- driving model assumptions used by the current engine
- solver and simulation settings
"""
function collect_run_parameters(intrinsic::IntrinsicConfig, job::JobConfig)
    col = get_column_parameters(condition=intrinsic.condition)
    a_vals, A_vals, C_vals = build_synapse_parameter_arrays(col)
    driving = get_default_driving_model()

    return Dict(
        "metadata" => Dict(
            "saved_at" => Dates.format(now(), "yyyy-mm-ddTHH:MM:SS"),
            "format" => "TOML",
            "note" => "Complete run snapshot for LaNMM Julia pipeline"
        ),
        "intrinsic_user" => to_dict(intrinsic),
        "column_parameters_effective" => col,
        "synapse_arrays_effective" => Dict(
            "a_vals" => a_vals,
            "A_vals" => A_vals,
            "C_vals" => C_vals
        ),
        "driving_model" => Dict(
            "e1" => merge(to_dict(driving.e1), Dict("baseline" => "swept via mu_p1_values")),
            "e2" => merge(to_dict(driving.e2), Dict("baseline" => "swept via mu_p2_values")),
            "pv" => to_dict(driving.pv),
            "seed" => driving.seed,
            "nonnegative_clipping" => driving.nonnegative_clipping
        ),
        "simulation" => Dict(
            "tmax_s" => job.tmax,
            "dt_s" => job.dt,
            "discard_s" => job.discard,
            "solver_primary" => "Tsit5",
            "solver_fallback" => "Rosenbrock23",
            "reltol" => 1e-8,
            "abstol" => 1e-8,
            "maxiters" => 10^8
        ),
        "analysis" => Dict(
            "alpha_band_hz" => collect(job.alpha_band),
            "gamma_band_hz" => collect(job.gamma_band),
            "power_method" => "bandpass_hilbert"
        ),
        "sweep" => Dict(
            "job_title" => job.job_title,
            "mu_p1_values" => job.mu_p1_values,
            "mu_p2_values" => job.mu_p2_values,
            "quiet_progress" => job.quiet_progress,
            "total_pairs" => length(job.mu_p1_values) * length(job.mu_p2_values)
        ),
        "nominal_run_for_diagnostics" => Dict(
            "mu_e1" => 270.0,
            "mu_e2" => 90.0,
            "xlims_plot" => [5, 10]
        )
    )
end

function save_parameters_snapshot(path::AbstractString, intrinsic::IntrinsicConfig, job::JobConfig)
    params = collect_run_parameters(intrinsic, job)
    open(path, "w") do io
        TOML.print(io, _to_toml_value(params); sorted=true)
    end
    return nothing
end

function _validate_intrinsic_config(intrinsic::IntrinsicConfig)
    intrinsic.fmax > 0 || error("intrinsic.fmax must be > 0")
    intrinsic.r_slope > 0 || error("intrinsic.r_slope must be > 0")
end

function _validate_job_config(job::JobConfig)
    !isempty(job.mu_p1_values) || error("job.mu_p1_values cannot be empty")
    !isempty(job.mu_p2_values) || error("job.mu_p2_values cannot be empty")

    job.dt > 0 || error("job.dt must be > 0")
    job.tmax > 0 || error("job.tmax must be > 0")
    job.discard >= 0 || error("job.discard must be >= 0")
    job.discard < job.tmax || error("job.discard must be < job.tmax")

    alpha = job.alpha_band
    gamma = job.gamma_band
    alpha[1] < alpha[2] || error("job.alpha_band must satisfy low < high")
    gamma[1] < gamma[2] || error("job.gamma_band must satisfy low < high")
end

function validate_configs(intrinsic::IntrinsicConfig, job::JobConfig)
    _validate_intrinsic_config(intrinsic)
    _validate_job_config(job)
    # Reuse condition whitelist from parameter layer
    get_column_parameters(condition=intrinsic.condition)
    return nothing
end

"""
Emit one progress update line for a sweep.
"""
function _print_progress(done::Int, total::Int, t_start::Float64; quiet::Bool=false)
    frac = done / total
    pct = 100 * frac
    elapsed = time() - t_start
    rate = done / max(elapsed, 1e-9)
    eta = (total - done) / max(rate, 1e-9)
    avg_per_sim = elapsed / max(done, 1)
    finish_clock = _fmt_finish_time(time() + eta)

    width = 30
    filled = clamp(round(Int, width * frac), 0, width)
    bar = repeat("=", filled) * repeat(".", width - filled)

    if quiet
        @printf(
            "Sweep progress [%s] %6.2f%% (%d/%d) | elapsed %s | ETA %s | avg %.3fs/sim | finish ~%s\n",
            bar, pct, done, total, _fmt_seconds(elapsed), _fmt_seconds(eta), avg_per_sim, finish_clock
        )
    else
        @printf(
            "\rSweep progress [%s] %6.2f%% (%d/%d) | elapsed %s | ETA %s | avg %.3fs/sim | finish ~%s",
            bar, pct, done, total, _fmt_seconds(elapsed), _fmt_seconds(eta), avg_per_sim, finish_clock
        )
        flush(stdout)
        if done == total
            println()
        end
    end
end

"""
Run a complete LaNMM sweep job.

Responsibilities:
- one nominal run + diagnostic plots
- full `(mu_p1, mu_p2)` sweep
- coupling/power/PEIX/frequency analyses
- figure and results serialization
"""
function run_sweep_job(intrinsic::IntrinsicConfig, job::JobConfig)
    validate_configs(intrinsic, job)

    job_title = job.job_title
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    output_dir = "$(job_title)_$(timestamp)"
    mkpath(output_dir)

    println("Saving parameters to $output_dir/parameters.txt")
    save_parameters_snapshot(joinpath(output_dir, "parameters.txt"), intrinsic, job)

    nominal = run_unified_simulation(
        intrinsic,
        job;
        mu_e1=270.0,
        mu_e2=90.0
    )
    plot_sim_results(nominal; save_dir=output_dir, xlims=(5, 10))

    m1_vals = job.mu_p1_values
    m2_vals = job.mu_p2_values
    total_pairs = length(m1_vals) * length(m2_vals)
    quiet_progress = job.quiet_progress
    update_every = quiet_progress ? max(1, total_pairs ÷ 10) : max(1, total_pairs ÷ 200)

    println("Running sweep for $(total_pairs) parameter pairs...")
    sweep_results = Dict{Tuple{Float64,Float64},Any}()
    done_pairs = 0
    t0 = time()

    for m1 in m1_vals, m2 in m2_vals
        sweep_results[(m1, m2)] = run_unified_simulation(intrinsic, job; mu_e1=m1, mu_e2=m2)
        done_pairs += 1
        if done_pairs == 1 || done_pairs % update_every == 0 || done_pairs == total_pairs
            _print_progress(done_pairs, total_pairs, t0; quiet=quiet_progress)
        end
    end

    println("Computing analysis products...")
    couplings = analyze_sweep_couplings(
        sweep_results;
        alpha_band=job.alpha_band,
        gamma_band=job.gamma_band
    )
    power_results = analyze_sweep_power(
        sweep_results;
        alpha_band=job.alpha_band,
        gamma_band=job.gamma_band,
        method=:bandpass_hilbert
    )
    peix_results = sweep_peix(
        sweep_results;
        r_slope=intrinsic.r_slope,
        v0_default=intrinsic.v0_default,
        v0_p2=intrinsic.v0_p2
    )
    freq_data = analyze_peak_frequencies(sweep_results; fs=1.0 / job.dt)

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
