using DifferentialEquations
using DataInterpolations
using DSP
using Random
using Statistics
if !isdefined(@__MODULE__, :IntrinsicConfig)
    include("LaNMM_Types.jl")
end
if !isdefined(@__MODULE__, :get_column_parameters)
    include("LaNMM_parameters.jl")
end

"""
Public API (this file)
----------------------
- `LaNMMParams`
- `build_intrinsic_params(; condition)`
- `run_unified_simulation(intrinsic, job; ...)`
- `run_unified_simulation(; ...)` (compatibility wrapper)

Internal helpers:
- `generate_am_signal`, `sigmoid_pop`, `lanmm_ode!`
"""

"""
Container for all parameters required by the ODE function.

The interpolation fields are parametric so the struct remains type-stable even if
the interpolation implementation changes.
"""
Base.@kwdef struct LaNMMParams{T1,T2,T3}
    v0_default::Float64
    v0_p2::Float64
    fmax::Float64
    r_slope::Float64
    a_vals::Vector{Float64}
    A_vals::Vector{Float64}
    C_vals::Vector{Float64}
    e1_interp::T1
    e2_interp::T2
    pv_interp::T3
end

"""
Build intrinsic synaptic arrays for a given condition.

Returns `(a_vals, A_vals, C_vals)` for synapses 1..14.
"""
function build_intrinsic_params(; condition::String="healthy")
    col = get_column_parameters(condition=condition)
    return build_synapse_parameter_arrays(col)
end

"""
Generate a simple amplitude-modulated input signal centered around `baseline`.
"""
function generate_am_signal(
    t_array::AbstractVector,
    carrier_freq,
    carrier_amplitude,
    baseline;
    envelope_mod_freq_hz=10.0,
    envelope_mod_depth=0.5
)
    envelope = carrier_amplitude .* (1.0 .+ envelope_mod_depth .* sin.(2π .* envelope_mod_freq_hz .* t_array))
    carrier = cos.(2π .* carrier_freq .* t_array)
    s = envelope .* carrier
    return (s .- mean(s)) .+ baseline
end

"""
Generate multiscale drive by combining slow AR(1)-like drift and fast low-pass noise.
"""
function generate_multiscale_noise(
    t_array::AbstractVector,
    base_mean::Float64;
    seed::Int=42,
    fs::Float64=1000.0,
    slow_std::Float64=400.0,
    slow_alpha::Float64=0.99,
    fast_std::Float64=5.0,
    fast_cutoff::Float64=100.0,
    floor_val::Float64=1e-4
)
    rng = MersenneTwister(seed)
    n = length(t_array)
    slow = zeros(Float64, n)
    slow[1] = base_mean

    for i in 2:n
        step = randn(rng) * slow_std
        slow[i] = slow_alpha * slow[i - 1] + (1.0 - slow_alpha) * (base_mean + step)
    end

    fast_white = randn(rng, n) .* fast_std
    filt = digitalfilter(Lowpass(fast_cutoff), Butterworth(4); fs=fs)
    fast = filtfilt(filt, fast_white)
    return max.(slow .+ fast, floor_val)
end

generate_drive_signal(t_array, mu::Float64, cfg::AMDrivingConfig; fs::Float64, seed::Int) = generate_am_signal(
    t_array,
    cfg.carrier_freq_hz,
    cfg.carrier_amplitude,
    mu;
    envelope_mod_freq_hz=cfg.envelope_mod_freq_hz,
    envelope_mod_depth=cfg.envelope_mod_depth
)

generate_drive_signal(t_array, mu::Float64, cfg::ConstantDrivingConfig; fs::Float64, seed::Int) = fill(mu, length(t_array))

generate_drive_signal(t_array, mu::Float64, cfg::MultiscaleDrivingConfig; fs::Float64, seed::Int) = generate_multiscale_noise(
    t_array,
    mu;
    seed=seed + cfg.seed_offset,
    fs=fs,
    slow_std=cfg.slow_std,
    slow_alpha=cfg.slow_alpha,
    fast_std=cfg.fast_std,
    fast_cutoff=cfg.fast_cutoff,
    floor_val=cfg.floor_val
)

@inline sigmoid_pop(v, v0, p::LaNMMParams) = p.fmax / (1.0 + exp(p.r_slope * (v0 - v)))

"""
LaNMM right-hand side with 14 second-order synapses (28 states).
"""
function lanmm_ode!(du, u, p::LaNMMParams, t)
    phi_e1 = max(0.0, p.e1_interp(t))
    phi_e2 = max(0.0, p.e2_interp(t))
    phi_pv = max(0.0, p.pv_interp(t))

    @inbounds begin
        vP1 = u[1] + u[3] + u[5] + u[21]
        vSS = u[7]
        vSST = u[9]
        vP2 = u[11] + u[13] + u[15] + u[23]
        vPV = u[17] + u[19] + u[25] + u[27]

        sP1 = sigmoid_pop(vP1, p.v0_default, p)
        sSS = sigmoid_pop(vSS, p.v0_default, p)
        sSST = sigmoid_pop(vSST, p.v0_default, p)
        sP2 = sigmoid_pop(vP2, p.v0_p2, p)
        sPV = sigmoid_pop(vPV, p.v0_default, p)

        presyn = (
            sSS, sSST, phi_e1, sP1, sP1, sP2, sPV,
            phi_e2, sP2, sPV, sP2, sP1, sP1, phi_pv
        )

        for s in 1:14
            idx_u = 2s - 1
            idx_z = 2s
            du[idx_u] = u[idx_z]
            du[idx_z] = p.a_vals[s] * p.A_vals[s] * (p.C_vals[s] * presyn[s]) -
                        2.0 * p.a_vals[s] * u[idx_z] -
                        (p.a_vals[s]^2) * u[idx_u]
        end
    end
    return nothing
end

"""
Run the unified simulation for one `(mu_e1, mu_e2)` pair.

Inputs:
- `intrinsic::IntrinsicConfig`: sigmoid and condition setup
- `job::JobConfig`: simulation timing config (`tmax`, `dt`, `discard`)
- keyword `mu_e1`, `mu_e2`, `mu_pv`: drive baselines
- keyword `driving::DrivingConfig`: driving waveform defaults

Output:
- `SimulationResult`
"""
function run_unified_simulation(
    intrinsic::IntrinsicConfig,
    job::JobConfig;
    mu_e1::Float64=270.0,
    mu_e2::Float64=90.0,
    mu_pv::Float64=0.0,
    driving::DrivingConfig=get_default_driving_model()
)
    # Build the simulation clock first; all generated drives are sampled on this grid.
    t_array = collect(0:job.dt:job.tmax)
    fs = 1.0 / job.dt
    a_vals, A_vals, C_vals = build_intrinsic_params(condition=intrinsic.condition)

    # Each population drive is generated independently from the selected mode
    # (constant / AM / multiscale), then optional clipping enforces nonnegative rates.
    e1_array = generate_drive_signal(t_array, mu_e1, driving.e1; fs=fs, seed=driving.seed)
    e2_array = generate_drive_signal(t_array, mu_e2, driving.e2; fs=fs, seed=driving.seed + 1)
    pv_array = generate_drive_signal(t_array, mu_pv, driving.pv; fs=fs, seed=driving.seed + 2)

    if driving.nonnegative_clipping
        e1_array = max.(0.0, e1_array)
        e2_array = max.(0.0, e2_array)
        pv_array = max.(0.0, pv_array)
    end

    p = LaNMMParams(
        v0_default=intrinsic.v0_default,
        v0_p2=intrinsic.v0_p2,
        fmax=intrinsic.fmax,
        r_slope=intrinsic.r_slope,
        a_vals=a_vals,
        A_vals=A_vals,
        C_vals=C_vals,
        e1_interp=LinearInterpolation(e1_array, t_array),
        e2_interp=LinearInterpolation(e2_array, t_array),
        pv_interp=LinearInterpolation(pv_array, t_array)
    )

    # 14 second-order synapses -> 28 ODE states.
    u0 = zeros(Float64, 28)
    prob = ODEProblem(lanmm_ode!, u0, (0.0, job.tmax), p)
    sol = solve(prob, Tsit5(), saveat=job.dt, reltol=1e-8, abstol=1e-8, maxiters=10^8)
    if !SciMLBase.successful_retcode(sol)
        # Fallback to a stiff-capable solver for hard parameter combinations.
        sol = solve(prob, Rosenbrock23(), saveat=job.dt, reltol=1e-8, abstol=1e-8, maxiters=10^8)
    end
    if !SciMLBase.successful_retcode(sol)
        error("ODE solve failed with retcode=$(sol.retcode). Try reducing dt/tmax or switching condition/drive settings.")
    end

    idx_start = searchsortedfirst(sol.t, job.discard)
    if idx_start > length(sol.t)
        error("Discard time ($(job.discard) s) is beyond the simulated interval ending at $(sol.t[end]) s.")
    end
    t_out = sol.t[idx_start:end]
    n = length(t_out)

    # Convert vector-of-state-vectors output to a dense matrix for fast slicing.
    u_mat = Matrix{Float64}(undef, n, 28)
    for (k, uk) in enumerate(sol.u[idx_start:end])
        @inbounds u_mat[k, :] = uk
    end

    # State mapping follows the Python model convention:
    # u(s) corresponds to first state of synapse s (odd positions in the 28-state vector).
    u_syn(s) = @view u_mat[:, 2s - 1]
    vP1 = u_syn(1) .+ u_syn(2) .+ u_syn(3) .+ u_syn(11)
    vP2 = u_syn(6) .+ u_syn(7) .+ u_syn(8) .+ u_syn(12)
    vPV = u_syn(9) .+ u_syn(10) .+ u_syn(13) .+ u_syn(14)

    # Re-evaluate drives exactly on simulated output times so plotting arrays
    # stay aligned even if the integrator terminates early.
    e1_out = [p.e1_interp(tt) for tt in t_out]
    e2_out = [p.e2_interp(tt) for tt in t_out]
    pv_out = [p.pv_interp(tt) for tt in t_out]

    return SimulationResult(
        t=t_out,
        vP1=collect(vP1),
        vP2=collect(vP2),
        vPV=collect(vPV),
        e1_array=e1_out,
        e2_array=e2_out,
        pv_array=pv_out,
        fs=fs,
        u1=collect(u_syn(1)),
        u2=collect(u_syn(2)),
        u3=collect(u_syn(3)),
        u4=collect(u_syn(4)),
        u5=collect(u_syn(5)),
        u6=collect(u_syn(6)),
        u7=collect(u_syn(7)),
        u8=collect(u_syn(8)),
        u9=collect(u_syn(9)),
        u10=collect(u_syn(10)),
        u11=collect(u_syn(11)),
        u12=collect(u_syn(12)),
        u13=collect(u_syn(13)),
        u14=collect(u_syn(14))
    )
end

"""
Compatibility wrapper using keyword-only configuration.
"""
function run_unified_simulation(;
    condition::String="healthy",
    mu_e1::Float64=270.0,
    mu_e2::Float64=90.0,
    mu_pv::Float64=0.0,
    v0_default::Float64=6.0,
    v0_p2::Float64=1.0,
    fmax::Float64=5.0,
    r_slope::Float64=0.56,
    tmax::Float64=4.0,
    dt::Float64=0.001,
    discard::Float64=1.0
)
    intrinsic = IntrinsicConfig(
        condition=condition,
        v0_default=v0_default,
        v0_p2=v0_p2,
        fmax=fmax,
        r_slope=r_slope
    )
    job = JobConfig(tmax=tmax, dt=dt, discard=discard)
    return run_unified_simulation(intrinsic, job; mu_e1=mu_e1, mu_e2=mu_e2, mu_pv=mu_pv)
end