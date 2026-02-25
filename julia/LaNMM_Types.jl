"""
Strongly-typed configuration and result objects for the Julia LaNMM pipeline.

Public API:
- `IntrinsicConfig`
- `JobConfig`
- `AMDrivingConfig`, `ConstantDrivingConfig`, `DrivingConfig`
- `SimulationResult`
- `SweepMetrics`
- `to_dict(...)` helpers for serialization/reporting
"""

Base.@kwdef struct IntrinsicConfig
    condition::String = "healthy"
    v0_default::Float64 = 6.0
    v0_p2::Float64 = 1.0
    fmax::Float64 = 5.0
    r_slope::Float64 = 0.56
end

Base.@kwdef struct JobConfig
    job_title::String = "julia_lanmm_sweep"
    mu_p1_values::Vector{Float64} = collect(50.0:20.0:400.0)
    mu_p2_values::Vector{Float64} = collect(50.0:20.0:400.0)
    tmax::Float64 = 10.0
    dt::Float64 = 0.001
    discard::Float64 = 1.0
    alpha_band::NTuple{2,Float64} = (8.0, 12.0)
    gamma_band::NTuple{2,Float64} = (30.0, 50.0)
    quiet_progress::Bool = false
end

Base.@kwdef struct AMDrivingConfig
    mode::Symbol = :am_simplified
    carrier_freq_hz::Float64 = 10.0
    carrier_amplitude::Float64 = 400.0
    envelope_mod_freq_hz::Float64 = 10.0
    envelope_mod_depth::Float64 = 0.5
end

Base.@kwdef struct MultiscaleDrivingConfig
    mode::Symbol = :multiscale
    seed_offset::Int = 0
    slow_std::Float64 = 400.0
    slow_alpha::Float64 = 0.99
    fast_std::Float64 = 5.0
    fast_cutoff::Float64 = 100.0
    floor_val::Float64 = 1e-4
end

Base.@kwdef struct ConstantDrivingConfig
    mode::Symbol = :constant
    baseline::Float64 = 0.0
end

Base.@kwdef struct DrivingConfig{T1,T2,T3}
    e1::T1 = AMDrivingConfig()
    e2::T2 = AMDrivingConfig(carrier_freq_hz=40.0)
    pv::T3 = ConstantDrivingConfig()
    seed::Int = 42
    nonnegative_clipping::Bool = true
end

"""
Typed simulation payload returned by `run_unified_simulation`.
"""
Base.@kwdef struct SimulationResult
    t::Vector{Float64}
    vP1::Vector{Float64}
    vP2::Vector{Float64}
    vPV::Vector{Float64}
    e1_array::Vector{Float64}
    e2_array::Vector{Float64}
    pv_array::Vector{Float64}
    fs::Float64
    u1::Vector{Float64}
    u2::Vector{Float64}
    u3::Vector{Float64}
    u4::Vector{Float64}
    u5::Vector{Float64}
    u6::Vector{Float64}
    u7::Vector{Float64}
    u8::Vector{Float64}
    u9::Vector{Float64}
    u10::Vector{Float64}
    u11::Vector{Float64}
    u12::Vector{Float64}
    u13::Vector{Float64}
    u14::Vector{Float64}
end

"""
Lightweight per-grid-point outputs for sweep runs.

Stores only scalar metrics so raw time-series can be discarded immediately.
"""
Base.@kwdef struct SweepMetrics
    r_s2e::Float64
    r_e2e::Float64
    p1_alpha::Float64
    p1_gamma::Float64
    p2_alpha::Float64
    p2_gamma::Float64
    peix_P1::Float64
    peix_P2::Float64
    alpha_peak_hz::Float64
    gamma_peak_hz::Float64
end

to_dict(x::IntrinsicConfig) = Dict(
    "condition" => x.condition,
    "v0_default" => x.v0_default,
    "v0_p2" => x.v0_p2,
    "fmax" => x.fmax,
    "r_slope" => x.r_slope
)

to_dict(x::JobConfig) = Dict(
    "job_title" => x.job_title,
    "mu_p1_values" => x.mu_p1_values,
    "mu_p2_values" => x.mu_p2_values,
    "tmax" => x.tmax,
    "dt" => x.dt,
    "discard" => x.discard,
    "alpha_band" => collect(x.alpha_band),
    "gamma_band" => collect(x.gamma_band),
    "quiet_progress" => x.quiet_progress
)

to_dict(x::AMDrivingConfig) = Dict(
    "mode" => string(x.mode),
    "carrier_freq_hz" => x.carrier_freq_hz,
    "carrier_amplitude" => x.carrier_amplitude,
    "envelope_mod_freq_hz" => x.envelope_mod_freq_hz,
    "envelope_mod_depth" => x.envelope_mod_depth
)

to_dict(x::MultiscaleDrivingConfig) = Dict(
    "mode" => string(x.mode),
    "seed_offset" => x.seed_offset,
    "slow_std" => x.slow_std,
    "slow_alpha" => x.slow_alpha,
    "fast_std" => x.fast_std,
    "fast_cutoff" => x.fast_cutoff,
    "floor_val" => x.floor_val
)

to_dict(x::ConstantDrivingConfig) = Dict(
    "mode" => string(x.mode),
    "baseline" => x.baseline
)

to_dict(x::DrivingConfig) = Dict(
    "e1" => to_dict(x.e1),
    "e2" => to_dict(x.e2),
    "pv" => to_dict(x.pv),
    "seed" => x.seed,
    "nonnegative_clipping" => x.nonnegative_clipping
)
